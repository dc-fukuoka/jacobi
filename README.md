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
only CPU[s]: 94.7290  
CPU + GPU[s]: 5.29180  
  
the following is an example of calculation result.
------------
m=1024 n=1024  
with GPU:  
~~~~
$ mpirun -x LD_LIBRARY_PATH -np 8 -npernode 4 -bind-to socket ./jacobi_acc_mpi.pgi
 iter_max:        10000
 m:         4096 n:         4096
 np:            8
 np_m:            1 np_n:            8
    0  0.499878
 1000  0.000357
 2000  0.000178
 3000  0.000118
 4000  0.000088
 5000  0.000071
 6000  0.000059
 7000  0.000050
 8000  0.000044
 9000  0.000039
10000  0.000035
Total CPU time:    5.29180
Device init   :    1.91826
Setup problem :   0.841110E-01
Data copyin   :   0.291271E-01
Computation   :    3.13129
Data copyout  :   0.277512E-01
Output        :   0.101262
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
