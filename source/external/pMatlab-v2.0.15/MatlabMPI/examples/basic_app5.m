%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic example of a typical application for which 
% MatlabMPI might be used.  It uses a leader that manages, 
% and workers that conduct the processing. This example 
% takes a matrix, breaks it up among processors, does a 
% computation, and gathers the results to display in a 
% figure window.
%
% This example builds on basic_app4.m. It keeps track of 
% each of the jobs that must be completed, and keeps track 
% of whether each of the worker processes are busy or not. 
% If there are more jobs to be completed and one of the worker 
% processes is not busy, then it sends the next sequential 
% row of the input matrix to that process. 
% It also asynchronously receives completed columns of the 
% output matrix in that the leader uses the MPI_Probe 
% function to determine whether a worker process has 
% sent a completed computation packet. The leader sends 
% a packet to be processed and then probes for a packet 
% that was sent to itself. If a packet has 
% been sent back, the leader processes the packet; if no 
% packet has been sent, then the leader determines whether 
% another work column should be sent to another worker. 
% If all of the workers are busy, the leader just waits 
% until any competed work has been returned.  
%
% This example can also plot the data on a graph figure. 
% To enable this feature, set the variable UseGraphics 
% equal to 1. 
%
% To run, start Matlab and type:
%
%   eval( MPI_Run('basic_app5',4,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('basic_app5',4,{'machine1' 'machine2' 'machine3' 'machine4'}) );
%
% Output will be piped into 4 files:
%
%   MatMPI/basic_app5.0.out
%   MatMPI/basic_app5.1.out
%   MatMPI/basic_app5.2.out
%   MatMPI/basic_app5.3.out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Albert Reuther
% MIT Lincoln Laboratory
% reuther@ll.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set whether a graph should be generated
UseGraphics = 1;

% Initialize MPI.
MPI_Init;

% Create communicator.
comm = MPI_COMM_WORLD;

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);

% Since the leader only manages, there must be at least 2 processes
if comm_size <= 1
    error('Cannot be run with only one process');
end

% Print rank.
disp(['my_rank: ',num2str(my_rank)]);

% Wait momentarily.
pause(1.0);

% Set who is the leader
leader = 0;
num_workers = comm_size - 1;

% Create base message tags.
coefs_tag = 10000;
input_tag = 20000;
output_tag = 30000;
end_tag = 999999;

% Set data sizes.
N1 = 1024;
N2 = 128;

% Leader.
if (my_rank == leader)
    % Create coefficient data - simple impluse.
    coefs = zeros(N1,1);
    coefs(1) = 1;
    
    % Create input data.
    input = ones(N1,N2);
    
    % Create output data array.
    output = zeros(N1,N2);
    
    if UseGraphics
        % Show output in a figure
        Hf_fig = imagesc(output);
        drawnow;
    end
    
    % Broadcast coefficients to everyone else.
    MPI_Bcast( leader, coefs_tag, comm, coefs );
    
    % flag for being done with all processing
    done = 0;
    
    % Instead of using for loops, use counters
    sendCounter = 1;
    recvCounter = 1;
    % Keep a vector of flags for whether each of the workers are busy ... 
    workers_busy = zeros(1, num_workers);
    % ... and whether each of the data items have been processed
    data_processed = zeros(1, N2);
    
    while ~done
        % Deal input data to everyone else (excluding self-leader).
        if sendCounter <= N2 
            if length(find(workers_busy)) < num_workers
                % Do not include leader in data dealing
                idle_workers = find(workers_busy == 0);
                dest = idle_workers(1);
                dest_tag = input_tag + sendCounter;
                dest_data = input(:,sendCounter);
                
                MPI_Send(dest,dest_tag,comm,dest_data);
                disp(['Sent data packet number ' num2str(sendCounter) ' to ' num2str(dest) ]);
                sendCounter = sendCounter + 1;
                workers_busy(dest) = 1;
            end
        end
        
        % Leader receives all the results.
        if length(find(data_processed)) < N2 
            
            [message_ranks, message_tags] = MPI_Probe( '*', '*', comm );
            % if message_ranks is not empty then receive the message
            if ~isempty(message_ranks)
                % Receive output.
                %disp(['Waiting on unit ' num2str(recvCounter)]);
                
                % Use this sort to try to service the received messages as close to actual order as possible
                [m_tags_sorted m_tags_sorted_idx] = sort(message_tags);
                dest = message_ranks(m_tags_sorted_idx(1));
                leader_tag = message_tags(m_tags_sorted_idx(1));
                mesg_num = leader_tag - output_tag;
                output(:, mesg_num) =  MPI_Recv( dest, leader_tag, comm);
                disp(['Received data packet number ' num2str(mesg_num) ' from ' num2str(dest)]);
                
                if UseGraphics
                    set(Hf_fig, 'CData', output);
                    drawnow;
                end
                data_processed(mesg_num) = 1;
                workers_busy(dest) = 0;
            else % is ~empty
                %disp(['Waiting on data packet ' num2str(recvCounter)]);
            end
            
        else    % recvCounter > N2
            done = 1;
        end
    end
    % When all of the work has been done, send out a finish broadcast so that 
    %   each of the worker processes exit properly
    MPI_Bcast(leader, end_tag, comm, done);
    
    disp('Sent ending broadcast');
    %end
    
    
    % Everyone but the leader receives the coefs.
else % (my_rank ~= leader)
    % Receive coefs.
    coefs = MPI_Recv( leader, coefs_tag, comm );
    disp('Received coefficients');
    done = 0;
    
    % Everyone but leader receives the input data and processes the results.
    while ~done
        [message_ranks, message_tags] = MPI_Probe( '*', '*', comm );
        % if message_ranks is not empty then receive the message
        if ~isempty(message_ranks)
            % Receive input.
            %disp(['Waiting on unit ' num2str(recvCounter)]);
            
            % if there are multiple messages waiting, make sure that 
            %  the end_tag is the last one read, since it will shut down
            %  the worker process - source should be leader (0)
            source = message_ranks(1);
            in_tag = message_tags(1);
            
            if in_tag ~= end_tag
                mesg_num = in_tag - input_tag;
                leader_tag = output_tag + mesg_num;
                
                % Receive input.
                i_input =  MPI_Recv(leader,in_tag,comm);
                disp(['Received data packet number ' num2str(mesg_num) ' from ' num2str(source)]);
                
                i_input = i_input + my_rank;
                
                % Do computation.
                i_output = fft(coefs) .* i_input;
                
                % Send results back to the leader.
                MPI_Send(leader,leader_tag,comm,i_output);
                
                disp(['Processed data unit ' num2str(mesg_num)]);
            else
                done = MPI_Recv(leader, end_tag, comm);
            end
            % If no message is waiting, microsleep to allow a message to be written
        else % is ~empty
            pause(0.01);
        end    
    end
end

% Finalize Matlab MPI.
disp('SUCCESS');
MPI_Finalize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2003 Massachusetts Institute of Technology
% 
% Permission is herby granted, without payment, to copy, modify, display
% and distribute this software and its documentation, if any, for any
% purpose, provided that the above copyright notices and the following
% three paragraphs appear in all copies of this software.  Use of this
% software constitutes acceptance of these terms and conditions.
%
% IN NO EVENT SHALL MIT BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF MIT HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 
% MIT SPECIFICALLY DISCLAIMS ANY EXPRESS OR IMPLIED WARRANTIES INCLUDING,
% BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.
%
% THIS SOFTWARE IS PROVIDED "AS IS," MIT HAS NO OBLIGATION TO PROVIDE
% MAINTENANCE, SUPPORT, UPDATE, ENHANCEMENTS, OR MODIFICATIONS.

