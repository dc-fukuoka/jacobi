jacobi - a benchmark by solving 2D laplace equation with jacobi iterative method.  
         GPU or Xeon Phi can be used.  
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
    # where $NP = np_m * np_n
    
    # if you do not have GPU nor Xeon Phi,
    $ mpirun -np $NP ./jacobi_mpi.intel.nooffload
    or
    $ mpirun -np $NP ./jacobi_mpi.pgi

    # to view the result
    $ gnuplot # splot "jacobi.dat" w l

performance comparison:
============
iter_max=10000  
m=4096, n=4096  
  
with Xeon E5-2680(8 cores/socket, 2 sockets) + Xeon Phi 5110P x2, 2 nodes (np_m=1, np_n=4):
------------
only CPU[s]: 66.7373  
CPU + MIC[s]: 29.1181  
  
with Xeon E5-2680 v4(14 cores/socket, 2 sockets) + NVIDIA Tesla P100 x4, 2 nodes (np_m=1, np_n=8):
------------
only CPU[s]: 76.6943  
CPU + GPU[s]: 5.29180  
  
the following is an example of calculation result.
------------
m=1024 n=1024  
with GPU:  
~~~~
$ mpirun -npernode 6 ./jacobi_acc_mpi.pgi
 iter_max:        10000
 m:         1024 n:         1024
 np:           12
 np_m:            1 np_n:           12
    0  0.499512
 1000  0.000349
 2000  0.000172
 3000  0.000113
 4000  0.000084
 5000  0.000067
 6000  0.000055
 7000  0.000047
 8000  0.000041
 9000  0.000036
10000  0.000032
Total CPU time:    4.18754
Device init   :    1.36887
Setup problem :   0.860450
Data copyin   :   0.410668E-01
Computation   :    1.87990
Data copyout  :   0.349710E-02
Output        :   0.337530E-01
~~~~
![Alt text](jacobi.gpu.gif?raw=true "calculated by GPU")

  
with Xeon Phi:  
~~~~
$ srun --mpi=pmi2 -N2 -n4 -c8 --cpu_bind=cores -t0:05:00 ./jacobi_mpi.intel
 iter_max:       10000
 m:        1024 n:        1024
 np:           4
 np_m:           1 np_n:           4
 
Number of threads = 236
    0  0.499512
 1000  0.000347
 2000  0.000171
 3000  0.000113
 4000  0.000084
 5000  0.000066
 6000  0.000055
 7000  0.000047
 8000  0.000041
 9000  0.000036
10000  0.000032
Total CPU time:    6.98289
Device init   :   0.400000E-05
Setup problem :   0.640400E-02
Data copyin   :   0.208880E-01
Computation   :    6.82102
Data copyout  :   0.292500E-02
Output        :   0.131650
~~~~
![Alt text](jacobi.mic.gif?raw=true "calculated by MIC")
