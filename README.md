jacobi - a benchmark by solving 2D laplace equation with jacobi iterative method.  
         intel Xeon Phi and GPUs can be used.  
============
how to run:  
  
    # for intel (intel compiler and intelmpi are required)  
    $ make mic  
    
    # for GPU (PGI compiler and its MPI are required)  
    $ make gpu  
    
    $ cat fort.18  
    10000       ! iter_max  
    4096 4096   ! m, n  
    1 2         ! np_m, np_n  
    $ vi fort.18 <-- adjust the paramters  
    
    $ mpirun -np $NP ./jacobi_mpi.intel   # for Xeon phi  
    or  
    $ mpirun -np $NP ./jacobi_acc_mpi.pgi # for GPU  
    $ gnuplot # splot "jacobi.dat" w l

performance comparison:
============
iter_max=10000  
m=4096, n=4096  
  
with Xeon E5-2680(8 cores/socket, 2 sockets) + Xeon Phi 5110P x2, 2 nodes:
------------
only CPU[s]: 66.7373  
CPU + MIC[s]: 29.1181  
  
with Xeon E5-2650(10 cores/socket, 2 sockets) + NVIDIA Tesla K80 x3, 2 nodes:
------------
only CPU[s]: 37.0406  
CPU + GPU[s]: 7.62425
