%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic example of a typical application for which 
% MatlabMPI might be used.  It uses a leader that manages, 
% and workers that conduct the processing. This example 
% takes a matrix, breaks it up among processors, does a 
% computation, and gathers the results to display in a 
% figure window.
%
% This example is different from basic_app, basic_app2, 
% and basic_app3 in that the leader uses the MPI_Probe 
% function to determine whether a worker process has 
% sent a completed computation packet. The leader sends 
% a packet to be processed and then probes for a packet 
% that was sent to itself. If a packet has 
% been sent back, the leader processes the packet; if no 
% packet has been sent, then the leader proceeds to 
% sending another packet to be processed. If no other 
% packets to be processed need to be sent, then it loops 
% waiting to receive a processed packet. 
%
% This example can also plot the data on a graph figure. 
% To enable this feature, set the variable UseGraphics 
% equal to 1. 
%
% To run, start Matlab and type:
%
%   eval( MPI_Run('basic_app2',4,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('basic_app2',4,{'machine1' 'machine2' 'machine3' 'machine4'}) );
%
% Output will be piped into 4 files:
%
%   MatMPI/basic_app2.0.out
%   MatMPI/basic_app2.1.out
%   MatMPI/basic_app2.2.out
%   MatMPI/basic_app2.3.out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Albert Reuther
% MIT Lincoln Laboratory
% reuther@ll.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set whether a graph should be generated
UseGraphics = 0;

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

% Create base message tags.
coefs_tag = 10000;
input_tag = 20000;
output_tag = 30000;

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
    
    while ~done
        % Deal input data to everyone else (excluding self-leader).
        if sendCounter <= N2 
            % Do not include leader in data dealing
            dest = mod((sendCounter - 1),(comm_size-1)) + 1;
            dest_tag = input_tag + sendCounter;
            dest_data = input(:,sendCounter);
            
            MPI_Send(dest,dest_tag,comm,dest_data);
            disp(['Sent data packet number ' num2str(sendCounter)]);
            sendCounter = sendCounter + 1;
        end
        
        % Leader receives all the results.
        if recvCounter <= N2 
            % Compute who sent this message.
            % Do not include leader in data dealing
            dest = mod((recvCounter - 1),(comm_size-1)) + 1;
            
            leader_tag = output_tag + recvCounter;
            
            [message_ranks, message_tags] = MPI_Probe( dest, leader_tag, comm );
            % if message_ranks is not empty then receive the message
            if ~isempty(message_ranks)
                % Receive output.
                %disp(['Waiting on unit ' num2str(recvCounter)]);
                output(:,recvCounter) =  MPI_Recv( dest, leader_tag, comm);
                disp(['Received data packet number ' num2str(recvCounter)]);
            
				if UseGraphics
					set(Hf_fig, 'CData', output);
					drawnow;
				end
                recvCounter = recvCounter + 1;
            else % is ~empty
                %disp(['Waiting on data packet ' num2str(recvCounter)]);
            end

        else    % recvCounter > N2
            done = 1;
        end
    end
%end


% Everyone but the leader receives the coefs.
else % (my_rank ~= leader)
    % Receive coefs.
    coefs = MPI_Recv( leader, coefs_tag, comm );
    disp('Received coefficients');
    
    % Everyone but leader receives the input data and processes the results.
    for i=1:N2
        % Compute who the destination is for this message.
        % Do not include leader in data dealing
        dest = mod((i - 1),(comm_size-1)) + 1;
        
        % Check if this destination is me.
        if (my_rank == dest)
            % Compute tags.
            dest_tag = input_tag + i;
            leader_tag = output_tag + i;
            
            % Receive input.
            i_input =  MPI_Recv(leader,dest_tag,comm);
            
            i_input = i_input + my_rank;
            
            % Do computation.
            i_output = fft(coefs) .* i_input;
            
            % Send results back to the leader.
            MPI_Send(leader,leader_tag,comm,i_output);
            
            disp(['Processed data unit ' num2str(i)]);
            pause(0.1);
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

