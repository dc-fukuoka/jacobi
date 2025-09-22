all : jacobi_acc_mpi.pgi

mic: jacobi_mpi.intel jacobi_mpi.intel.nooffload
gpu: jacobi_mpi.pgi jacobi_acc_mpi.pgi

jacobi_mpi.intel: jacobi.f90
	mpiifort -D_MPI -DOVERLAP -o $@ -fpp -march=core-avx2 -O3 -fopenmp $<

jacobi_mpi.intel.nooffload: jacobi.f90
	mpiifort -D_MPI -DOVERLAP -o $@ -fpp -march=core-avx2 -O3 -fopenmp -qno-offload $<

jacobi_mpi.pgi: jacobi.f90
	mpif90 -D_MPI -DOVERLAP -o $@ -Mpreprocess -O4 -Mvect=simd:256 -Mfma -mp $<

jacobi_acc_mpi.pgi: jacobi.f90
	mpif90 -D_MPI -DOVERLAP -o $@ -Mpreprocess -O4 -Mvect=simd:256 -Mfma -acc -Minfo=accel $<

clean:
	rm -f jacobi_mpi.intel jacobi_mpi.intel.nooffload jacobi_mpi.pgi jacobi_acc_mpi.pgi jacobi.dat fort.999 *.mod *.modmic

clean-mic:
	rm -f jacobi_mpi.intel jacobi_mpi.intel.nooffload jacobi.dat fort.999 *.mod *.modmic

clean-gpu:
	rm -f jacobi_mpi.pgi jacobi_acc_mpi.pgi jacobi.dat fort.999 *.mod *.modmic
