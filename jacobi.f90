module params
#ifdef __INTEL_OFFLOAD
use mic_lib
#endif
  implicit none
#ifdef _MPI
  include "mpif.h"
#endif

  integer,parameter :: d_=8 ! floating point precision
  integer :: m,n ! matrix size

#ifdef _MPI
  integer :: iam, np, ierr
  integer :: ireq,iprov
  integer :: np_m, np_n
  integer :: up,down,left,right,coords(2)
  integer :: idiv(2)
  logical :: periods(2)
  integer :: comm_cart, iam_cart
  integer ::istat(mpi_status_size)
  integer :: istats_exchange(mpi_status_size,8)
  integer :: ireqs_exchange(8)
  integer,allocatable :: ireqs_dist(:),istats_dist(:,:)
  integer :: type_a, type_a_l_rest_n, type_a_l_rest_m, type_a_l_norest, type_updown
  real(kind=d_),dimension(:,:),allocatable,target :: a_l,anew_l
  real(kind=d_),dimension(:,:),allocatable        :: buf_m,buf_n
  integer :: m_l, n_l
  integer :: m_l_eq, n_l_eq
  integer :: mod_m, mod_n
  integer :: tmp_m, tmp_n
  logical :: if_update_left,if_update_right,if_update_up,if_update_down

#ifdef __INTEL_OFFLOAD
  real(kind=d_),dimension(:),pointer :: p_left, p_right, p_up, p_down
  integer :: nummics, mymic
  integer :: sig
  !dir$ attributes offload:mic :: m_l, n_l
  !dir$ attributes offload:mic :: p_left, p_right, p_up, p_down
  !dir$ attributes align:64    :: a_l,anew_l
  !dir$ attributes offload:mic :: a_l,anew_l
#endif

  integer,dimension(:),allocatable :: index_m, index_n
  real(kind=d_) :: gerror
#ifdef _OPENACC
  integer:: numdevices, mydevice
#endif
#endif
  integer :: i,j,k
  real(kind=d_),dimension(:,:),allocatable :: a,anew
  real(kind=d_),parameter                  :: tol=0.0_d_ ! error tolerance

  integer           :: iter_max ! maxmum iteration count
  real(kind=d_)     :: error
  integer           :: iter
  integer,parameter :: dskip=32
  integer           :: pr_iter
  real(kind=d_)     :: time(7)
  integer(8)        :: clock(7),clock_rate,max_count

#if defined(__INTEL_OFFLOAD) && !defined(_MPI)
  !dir$ attributes align:64    :: a,anew 
  !dir$ attributes offload:mic :: a,anew
  !dir$ attributes offload:mic :: i,j,m,n
#endif
  !dir$ attributes offload:mic :: error,i,j
end module params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program jacobi
  use params
  implicit none

  call initialize()

#ifdef _MPI
  call setup_mpi() 
#endif

  call initialize_device()

  call setup_problem()

#ifdef _MPI
  call distribute_data()
#endif

  call copyin()

  call jacobi_laplace2d()

  call copyout()

#ifdef _MPI
  call gather_data()
#endif

  call write_result()

  call diagnostics()

  call finalize()

  stop
end program jacobi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize()
  use params
  implicit none

#ifdef _MPI  
  ireq=mpi_thread_funneled
  call mpi_init_thread(ireq,iprov,ierr)
  if (iprov<ireq) then
     write(6,*) "mpi_thread_funneled is required."
     call mpi_finalize(ierr)
     stop
  end if

  call mpi_comm_rank(mpi_comm_world,iam,ierr)
  call mpi_comm_size(mpi_comm_world,np, ierr)
#endif

  read(18,*) iter_max
  read(18,*) m,n
  pr_iter=max(iter_max/10,1)

#ifdef _MPI
  read(18,*) np_m, np_n

  if (np.ne.np_m*np_n) then
     write(6,*) "np must equal to np_m*np_n=",np_m*np_n
     call mpi_finalize(ierr)
     stop
  end if

  if (iam.eq.0) then
#endif
     write(6,*) "iter_max:",iter_max
     write(6,*) "m:",m,"n:",n
#ifdef _MPI
     write(6,*) "np:", np
     write(6,*) "np_m:",np_m,"np_n:",np_n
  end if
#endif

  return
end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _MPI
subroutine setup_mpi()
  use params
  implicit none

  allocate(index_m(0:np-1), index_n(0:np-1))
  idiv(1) = np_m
  idiv(2) = np_n
  periods(:) = .false.

  call mpi_cart_create(mpi_comm_world,2,idiv,periods, &
       .false.,comm_cart,ierr)
  call mpi_comm_rank(comm_cart,iam_cart,ierr)
  call mpi_cart_coords(comm_cart,iam_cart,2,coords,ierr)

  mod_m = mod(m,np_m)
  mod_n = mod(n,np_n)
  if (coords(1).eq.0) then
     m_l = m / np_m + mod_m
  else
     m_l = m / np_m
  end if

  if (coords(2).eq.0) then
     n_l = n / np_n + mod_n
  else
     n_l = n / np_n
  end if

  m_l_eq = m / np_m
  n_l_eq = n / np_n

  allocate(a_l(0:m_l+1,0:n_l+1), anew_l(m_l,n_l))
  allocate(buf_m(1:m_l,1:2), buf_n(1:n_l,1:2))
  a_l(:,:) = 0.0d0
  anew_l(:,:) = 0.0d0
  buf_m(:,:) = 0.0d0
  buf_n(:,:) = 0.0d0

  ! create a datatype for halo data exchange
  call mpi_type_vector(n_l, 1, m_l+2, mpi_real8, type_updown, ierr)
  call mpi_type_commit(type_updown,ierr)

  call mpi_cart_shift(comm_cart,1,1,left,right,ierr)
  call mpi_cart_shift(comm_cart,0,1,up,  down, ierr)

  ! for distribution
  ! for a
  call mpi_type_vector(n_l+2, m_l+2, m+2, mpi_real8, type_a, ierr)
  call mpi_type_commit(type_a,  ierr)

  ! for process 0
  ! for a_l with the remainder of mod_n
  call mpi_type_vector(n_l_eq+mod_n+2, m_l_eq+2, m+2, mpi_real8, type_a_l_rest_n, ierr)
  call mpi_type_commit(type_a_l_rest_n,  ierr)

  ! for a_l with the remainder of mod_m
  call mpi_type_vector(n_l_eq+2, m_l_eq+mod_m+2, m+2, mpi_real8, type_a_l_rest_m, ierr)
  call mpi_type_commit(type_a_l_rest_m,  ierr)

  ! for a_l without any remainder
  call mpi_type_vector(n_l_eq+2, m_l_eq+2,       m+2, mpi_real8, type_a_l_norest, ierr)
  call mpi_type_commit(type_a_l_norest,  ierr)

!  calculate indexes
  k=0
  tmp_m=0
  tmp_n=0
  do j = 1, np_m
     do i = 1, np_n
        if (j.eq.2) tmp_m = mod_m
        if (i.eq.2) tmp_n = mod_n
        index_m(k) = 1 + (m/np_m)*(j-1) + tmp_m
        index_n(k) = 1 + (n/np_n)*(i-1) + tmp_n
        k = k + 1
     end do
     tmp_n = 0
  end do

! if the process is on the border of the process grid, there is no need to update halo data for the direction.
  if_update_left  = .true.
  if_update_right = .true.
  if_update_up    = .true.
  if_update_down  = .true.
  if (coords(2).eq.0)      if_update_left  = .false.
  if (coords(2).eq.np_n-1) if_update_right = .false.
  if (coords(1).eq.0)      if_update_up    = .false.
  if (coords(1).eq.np_m-1) if_update_down  = .false.

  return
end subroutine setup_mpi
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_device()
#ifdef _OPENACC
  use openacc
#endif
#ifdef _OPENMP
  !$ use omp_lib
#endif
  use params
  implicit none

#ifdef _MPI
  if (iam.eq.0) then
#endif 
#ifdef _OPENMP
#ifdef __INTEL_OFFLOAD
  !dir$ offload target(mic)
#endif
  !$omp parallel
  !$omp master
  !$ write(*,'("Number of threads = ",i3)') omp_get_num_threads()
  !$omp end master
  !$omp end parallel
#endif
#ifdef _MPI
  end if
#endif

  ! =================
  ! Initialyze device
  ! =================
  !call cpu_time(time(1))
  call system_clock(clock(1),clock_rate,max_count)
#if !defined(GCC) && defined(_OPENACC)
  !$acc init device_type(nvidia)
#endif
#ifdef GCC
  call acc_init(acc_device_nvidia)
#endif

#if defined(_MPI) && defined(_OPENACC)
  numdevices = acc_get_num_devices(acc_device_nvidia)
  mydevice = mod(iam, numdevices)
  call acc_set_device_num(mydevice,acc_device_nvidia)
#endif

#if defined(_MPI) && defined(__INTEL_OFFLOAD)
  nummics = offload_number_of_devices()
  mymic = mod(iam,nummics)
#endif

  return
end subroutine initialize_device
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setup_problem()
  use params
  implicit none

  ! =============
  ! Setup problem
  ! =============
  !call cpu_time(time(2))
  call system_clock(clock(2),clock_rate,max_count)

#ifdef _MPI
  if (iam.eq.0) then
#endif
     allocate(a(0:m+1,0:n+1),anew(m,n))
     ! Initial and boundary condition
     a(:,:)=0.0_d_

     do i=0,m+1
        a(i,0)=real(m+1-i,d_)/(m+1)
     end do
     do j=0,n+1
        a(0,j)=real(n+1-j,d_)/(n+1)
     end do
#ifdef _MPI
  end if
#endif

  return
end subroutine setup_problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _MPI
subroutine distribute_data
  use params
  implicit none
  ! distribute the data 
  if (iam.eq.0) then
     allocate(ireqs_dist(1:np-1))
     allocate(istats_dist(mpi_status_size,1:np-1))
     k=0
     do j = 0, np_m-1
        do i = 0, np_n-1
           if (i.eq.0.and.j.eq.0) then 
              a_l(0:m_l+1,0:n_l+1) = a(0:m_l+1,0:n_l+1)
           else if (i.eq.0.and.j.ne.0) then
             call mpi_isend(a(index_m(k)-1,index_n(k)-1), 1, type_a_l_rest_n, k, &
                  k, mpi_comm_world, ireqs_dist(k), ierr)
           else if (i.ne.0.and.j.eq.0) then
             call mpi_isend(a(index_m(k)-1,index_n(k)-1), 1, type_a_l_rest_m, k, &
                  k, mpi_comm_world, ireqs_dist(k), ierr)
           else if (i.ne.0.and.j.ne.0) then
             call mpi_isend(a(index_m(k)-1,index_n(k)-1), 1, type_a_l_norest, k, &
                  k, mpi_comm_world, ireqs_dist(k), ierr)
           end if
           k = k + 1
        end do
     end do
     call mpi_waitall(np-1,ireqs_dist,istats_dist,ierr)
     deallocate(ireqs_dist,istats_dist)
  else
     call mpi_recv(a_l(0,0), (m_l+2)*(n_l+2), mpi_real8, 0, &
          iam, mpi_comm_world,istat, ierr)
  end if
  
  return
end subroutine distribute_data
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine copyin()
  use params
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef _OPENACC
  use openacc
#endif
  implicit none

  ! ===========
  ! Data copyin
  ! ===========
  !call cpu_time(time(3))
  call system_clock(clock(3),clock_rate,max_count)
#if defined(_MPI) && defined(_OPENACC)
  !$acc enter data copyin(a_l(:,:)) create(anew_l(:,:))
#elif !defined(_MPI) && defined(_OPENACC)
  !$acc enter data copyin(a(:,:)) create(anew(:,:))
#endif

#if defined(_MPI) && defined(__INTEL_OFFLOAD)
  ! workaround
  p_left  => a_l(1:m_l,    0)
  p_right => a_l(1:m_l,n_l+1)
  p_up    => a_l(0,    1:n_l)
  p_down  => a_l(m_l+1,1:n_l)
  !dir$ offload_transfer target(mic:mymic)     in(a_l(:,:)    : alloc_if(.true.) free_if(.false.)) signal(sig)
  !dir$ offload_transfer target(mic:mymic) nocopy(anew_l(:,:) : alloc_if(.true.) free_if(.false.)) signal(sig)
  !dir$ offload_transfer target(mic:mymic) nocopy(p_left      : alloc_if(.true.) free_if(.false.)) signal(sig)
  !dir$ offload_transfer target(mic:mymic) nocopy(p_right     : alloc_if(.true.) free_if(.false.)) signal(sig)
  !dir$ offload_transfer target(mic:mymic) nocopy(p_up        : alloc_if(.true.) free_if(.false.)) signal(sig)
  !dir$ offload_transfer target(mic:mymic) nocopy(p_down      : alloc_if(.true.) free_if(.false.)) signal(sig)
  !dir$ offload_wait target(mic:mymic) wait(sig)
#endif

  return
end subroutine copyin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine jacobi_laplace2d()
  use params
  implicit none

  ! ===========
  ! Computation
  ! ===========
  !call cpu_time(time(4))
  call system_clock(clock(4),clock_rate,max_count)
  ! Jacobi method

#ifdef _MPI
  call mpi_barrier(mpi_comm_world,ierr)
#endif

  do iter=0,iter_max
     error=0.0_d_

#ifdef _MPI
! halo data exchange
     ! update send address on the host
     !$acc update if(if_update_left)  host(a_l(1:m_l,    1))
     !$acc update if(if_update_right) host(a_l(1:m_l,  n_l))
     !$acc update if(if_update_up)    host(a_l(1,    1:n_l))
     !$acc update if(if_update_down)  host(a_l(m_l,  1:n_l))

#ifdef __INTEL_OFFLOAD
     ! workaround
     if (if_update_left)  p_left  => a_l(1:m_l,  1)
     if (if_update_right) p_right => a_l(1:m_l,n_l)
     if (if_update_up)    p_up    => a_l(1,  1:n_l)
     if (if_update_down)  p_down  => a_l(m_l,1:n_l)
     if (if_update_left) then
        !dir$ offload_transfer target(mic:mymic) out(p_left  : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     if (if_update_right) then
        !dir$ offload_transfer target(mic:mymic) out(p_right : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     if (if_update_up) then
        !dir$ offload_transfer target(mic:mymic) out(p_up    : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     if (if_update_down) then
        !dir$ offload_transfer target(mic:mymic) out(p_down  : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     !dir$ offload_wait target(mic:mymic) wait(sig)
#endif

#ifdef OVERLAP_OMP
     !$omp parallel default(shared)
     !$omp master
#endif

! left, right
     call mpi_irecv(buf_m(1,1), m_l, mpi_real8,     right, 0, &
          comm_cart, ireqs_exchange(1), ierr)
     call mpi_irecv(buf_m(1,2), m_l, mpi_real8,     left,  1, &
          comm_cart, ireqs_exchange(2), ierr)

     call mpi_isend(a_l(1,1),   m_l, mpi_real8,     left,  0, &
          comm_cart, ireqs_exchange(3), ierr)
     call mpi_isend(a_l(1,n_l), m_l, mpi_real8,     right, 1, &
          comm_cart, ireqs_exchange(4), ierr)
! up, down
     call mpi_irecv(buf_n(1,1), n_l,   mpi_real8,   down,  2, &
          comm_cart, ireqs_exchange(5), ierr)
     call mpi_irecv(buf_n(1,2), n_l,   mpi_real8,    up,   3, &
          comm_cart, ireqs_exchange(6), ierr)

     call mpi_isend(a_l(1,1),       1, type_updown,  up,   2, &
          comm_cart, ireqs_exchange(7), ierr)
     call mpi_isend(a_l(m_l,1),     1, type_updown,  down, 3, &
          comm_cart, ireqs_exchange(8), ierr)

#ifndef OVERLAP_MPI
     ! update halos
     call mpi_waitall(8, ireqs_exchange, istats_exchange, ierr)
#endif
#endif

#if defined(__INTEL_OFFLOAD) && defined(_MPI) && !defined(OVERLAP_OMP)
     !dir$ offload target(mic:mymic) &
     !dir$ nocopy(anew_l(:,:) : alloc_if(.false.) free_if(.false.)) &
     !dir$ nocopy(a_l(:,:)    : alloc_if(.false.) free_if(.false.)) &
     !dir$ inout(error)
#elif defined(__INTEL_OFFLOAD) && !defined(_MPI)
     !dir$ offload target(mic)
#endif
#if defined(_OPENMP) && !defined(_OPENACC) && !defined(OVERLAP_OMP)
     !$omp parallel default(shared) private(i,j)
     !$omp do reduction(max:error)
#elif defined(OVERLAP_OMP)
     !$omp end master
     !$omp do reduction(max:error) private(i,j)
#endif
#ifdef _MPI
     !$acc parallel private(i,j) present(anew_l(:,:),a_l(:,:))
     !$acc loop collapse(2) reduction(max:error)
     do j=1,n_l
        !dir$ simd
        do i=1,m_l
           anew_l(i,j) = 0.25_d_*(a_l(i+1,j)+a_l(i-1,j) &
                +a_l(i,j-1)+a_l(i,j+1))
           error = max(error,abs(anew_l(i,j)-a_l(i,j)))
        end do
     end do
#else
     !$acc parallel private(i,j) present(anew(:,:),a(:,:))
     !$acc loop collapse(2) reduction(max:error)
     do j=1,n
        !dir$ simd
        do i=1,m
           anew(i,j) = 0.25_d_*(a(i+1,j)+a(i-1,j) &
                +a(i,j-1)+a(i,j+1))
           error = max(error,abs(anew(i,j)-a(i,j)))
        end do
     end do
#endif
#if defined(_OPENMP) && !defined(_OPENACC)
     !$omp end do nowait
     !$omp do
#endif
#ifdef _MPI
     !$acc loop collapse(2)
     do j=1,n_l
        !dir$ simd
        do i=1,m_l
           a_l(i,j) = anew_l(i,j)
        end do
     end do
#else
     !$acc loop collapse(2)
     do j=1,n
        !!dir$ simd
        do i=1,m
           a(i,j) = anew(i,j)
        end do
     end do     
#endif
     !$acc end parallel
#if defined(_OPENMP) && !defined(_OPENACC)
     !$omp end parallel
#endif

#ifdef _MPI
#ifdef OVERLAP_MPI
     ! update halo
     call mpi_waitall(8, ireqs_exchange, istats_exchange, ierr)
#endif

     call mpi_allreduce(error, gerror, 1, mpi_real8, mpi_max, &
          mpi_comm_world, ierr)
     error = gerror

     if (if_update_left)  a_l(1:m_l, 0    ) = buf_m(1:m_l,2)
     if (if_update_right) a_l(1:m_l, n_l+1) = buf_m(1:m_l,1)
     if (if_update_up)    a_l(0,     1:n_l) = buf_n(1:n_l,2)
     if (if_update_down)  a_l(m_l+1, 1:n_l) = buf_n(1:n_l,1)
     !$acc update if(if_update_left)  device(a_l(1:m_l,    0))
     !$acc update if(if_update_right) device(a_l(1:m_l,n_l+1))
     !$acc update if(if_update_up)    device(a_l(0,    1:n_l))
     !$acc update if(if_update_down)  device(a_l(m_l+1,1:n_l))

#ifdef __INTEL_OFFLOAD
     ! workaround
     if (if_update_left)  p_left  => a_l(1:m_l,    0)
     if (if_update_right) p_right => a_l(1:m_l,n_l+1)
     if (if_update_up)    p_up    => a_l(0,    1:n_l)
     if (if_update_down)  p_down  => a_l(m_l+1,1:n_l)
     if (if_update_left) then
        !dir$ offload_transfer target(mic:mymic) in(p_left  : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     if (if_update_right) then
        !dir$ offload_transfer target(mic:mymic) in(p_right : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     if (if_update_up) then
        !dir$ offload_transfer target(mic:mymic) in(p_up    : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     if (if_update_down) then
        !dir$ offload_transfer target(mic:mymic) in(p_down  : alloc_if(.false.) free_if(.false.)) signal(sig)
     end if
     !dir$ offload_wait target(mic:mymic) wait(sig)
#endif

     if (iam.eq.0) then
#endif
        if(error<tol) then
           write(*,'(i5,f10.6)') iter,error
           exit
        end if
        if(mod(iter,pr_iter)==0) write(*,'(i5,f10.6)') iter,error
#ifdef _MPI
     end if
#endif
  end do

#ifdef _MPI
  call mpi_barrier(mpi_comm_world,ierr)
#endif

  return
end subroutine jacobi_laplace2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine copyout()
  use params
  implicit none

  ! ============
  ! Data copyout
  ! ============
  !call cpu_time(time(5))
  call system_clock(clock(5),clock_rate,max_count)
#ifdef _MPI
  !$acc exit data copyout(a_l(:,:)) &
  !$acc delete(anew_l(:,:))
#else
  !$acc exit data copyout(a(:,:)) &
  !$acc delete(anew(:,:))
#endif

#if defined(_MPI) && defined(__INTEL_OFFLOAD)
  !dir$ offload_transfer target(mic:mymic) out(a_l(:,:)       : alloc_if(.false.) free_if(.true.)) signal(sig)
  !dir$ offload_transfer target(mic:mymic) nocopy(anew_l(:,:) : alloc_if(.false.) free_if(.true.)) signal(sig)
  !dir$ offload_wait target(mic:mymic) wait(sig)
  nullify(p_left,p_right,p_up,p_down)
#endif

  return
end subroutine copyout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _MPI
subroutine gather_data()
  use params
  implicit none

! data gather

! need to update the corner data before gather
! left, right
     call mpi_irecv(a_l(0,0),         1, mpi_real8, left,  0, &
          comm_cart, ireqs_exchange(1), ierr)
     call mpi_irecv(a_l(0,n_l+1),     1, mpi_real8, right, 1, &
          comm_cart, ireqs_exchange(2), ierr)
     call mpi_irecv(a_l(m_l+1,0),     1, mpi_real8, left,  2, &
          comm_cart, ireqs_exchange(3), ierr)
     call mpi_irecv(a_l(m_l+1,n_l+1), 1, mpi_real8, right, 3, &
          comm_cart, ireqs_exchange(4), ierr)

     call mpi_isend(a_l(0,n_l),     1, mpi_real8,   right, 0, &
          comm_cart, ireqs_exchange(5), ierr)
     call mpi_isend(a_l(0,1),       1, mpi_real8,   left,  1, &
          comm_cart, ireqs_exchange(6), ierr)
     call mpi_isend(a_l(m_l+1,n_l), 1, mpi_real8,   right, 2, &
          comm_cart, ireqs_exchange(7), ierr)
     call mpi_isend(a_l(m_l+1,1),   1, mpi_real8,   left,  3, &
          comm_cart, ireqs_exchange(8), ierr)
     call mpi_waitall(8, ireqs_exchange, istats_exchange, ierr)

  if (iam.eq.0) then
! debug
!     a(:,:) = -1.0d0

     k=0
     do j = 0, np_m-1
        do i = 0, np_n-1
           if (i.eq.0.and.j.eq.0) then
              a(0:m_l+1,0:n_l+1) = a_l(0:m_l+1,0:n_l+1)
           else if (i.eq.0.and.j.ne.0) then
             call mpi_recv(a(index_m(k)-1,index_n(k)-1), 1, type_a_l_rest_n, k, &
                  k, mpi_comm_world, istat, ierr)
           else if (i.ne.0.and.j.eq.0) then
             call mpi_recv(a(index_m(k)-1,index_n(k)-1), 1, type_a_l_rest_m, k, &
                  k, mpi_comm_world, istat, ierr)
           else if (i.ne.0.and.j.ne.0) then
             call mpi_recv(a(index_m(k)-1,index_n(k)-1), 1, type_a_l_norest, k, &
                  k, mpi_comm_world, istat, ierr)
           end if
           k = k + 1
        end do
     end do

  else
     call mpi_send(a_l(0,0), (m_l+2)*(n_l+2), mpi_real8, 0, iam, &
          mpi_comm_world, ierr)
  end if

  return
end subroutine gather_data
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_result()
  use params
  implicit none

  ! ======
  ! Output
  ! ======
  !call cpu_time(time(6))
  call system_clock(clock(6),clock_rate,max_count)

#ifdef _MPI
  if (iam.eq.0) then
#endif
     write(999) a
     open(10,file='jacobi.dat')
     do j=0,n+1,dskip
        do i=0,m+1,dskip
           write(10,*) i,j,a(i,j)
        end do
        if(mod(m+1,dskip)/=0) write(10,*) m+1,j,a(m+1,j)
        write(10,*)
     end do
     if(mod(n+1,dskip)/=0) then
        do i=0,m+1,dskip
           write(10,*) i,n+1,a(i,n+1)
        end do
        if(mod(m+1,dskip)/=0) write(10,*) m+1,n+1,a(m+1,n+1)
     end if
     close(10)
#ifdef _MPI
  end if
#endif

  !call cpu_time(time(7))
  call system_clock(clock(7),clock_rate,max_count)

  return
end subroutine write_result
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagnostics()
  use params
  implicit none

  ! ===========
  ! Diagnostics
  ! ===========
#ifdef _MPI
  if (iam.eq.0) then
#endif
     do i=1,6
        if(clock(i+1)>=clock(i)) then
           clock(i)=clock(i+1)-clock(i)
        else
           clock(i)=(max_count-clock(i))+clock(i+1)
        end if
        time(i)=real(clock(i),8)/clock_rate
     end do
     time(7)=sum(time(1:6))
     write(*,'("Total CPU time: ",g14.6)') time(7)
     write(*,'("Device init   : ",g14.6)') time(1)
     write(*,'("Setup problem : ",g14.6)') time(2)
     write(*,'("Data copyin   : ",g14.6)') time(3)
     write(*,'("Computation   : ",g14.6)') time(4)
     write(*,'("Data copyout  : ",g14.6)') time(5)
     write(*,'("Output        : ",g14.6)') time(6)
#ifdef _MPI
  end if
#endif

  return
end subroutine diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finalize()
  use params
  implicit none

#ifdef _MPI
  if (iam.eq.0) then
#endif
     deallocate(a,anew)
#ifdef _MPI
  end if

  deallocate(a_l,anew_l)
  deallocate(buf_m,buf_n)
  deallocate(index_m,index_n)

  call mpi_type_free(type_updown,       ierr)
  call mpi_type_free(type_a,            ierr)
  call mpi_type_free(type_a_l_rest_n,   ierr)
  call mpi_type_free(type_a_l_rest_m,   ierr)
  call mpi_type_free(type_a_l_norest,   ierr)
  call mpi_comm_free(comm_cart,         ierr)
  call mpi_finalize(ierr)
#endif

  return
end subroutine finalize
