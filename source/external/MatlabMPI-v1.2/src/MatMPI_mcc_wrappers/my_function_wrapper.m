function dummy=my_function(my_cpu)
global MPI_COMM_WORLD;
load 'MatMPI/MPI_COMM_WORLD.mat';
MPI_COMM_WORLD.rank = my_cpu;
