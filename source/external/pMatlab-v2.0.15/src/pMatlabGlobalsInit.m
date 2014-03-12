% Initialize globals necessary for pMatlab.
global pMATLAB
MPI_COMM_WORLD.rank = 0;
MPI_COMM_WORLD.size = 1;
MPI_COMM_WORLD.save_message_flag = 0;
MPI_COMM_WORLD.group = (1:MPI_COMM_WORLD.size)-1;
MPI_COMM_WORLD.machine_id = zeros(1,MPI_COMM_WORLD.size);
MPI_COMM_WORLD.host_rank = 0;
pMATLAB.comm = MPI_COMM_WORLD; %MPI communicator
pMATLAB.comm_size = MPI_COMM_WORLD.size;
pMATLAB.my_rank = MPI_COMM_WORLD.rank;
pMATLAB.leader=0; %set the leader
pMATLAB.host = MPI_COMM_WORLD.rank;
pMATLAB.pList = [0:pMATLAB.comm_size-1]; %set processor list
pMATLAB.tag_num = 0;
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

clear('MPI_COMM_WORLD');
