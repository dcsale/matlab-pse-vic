function sync_start(comm,starter_rank,delay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function tries to get all processes in a
% MatlabMPI communicatior to wait and start at the same.
% Assumes clocks across the machines are synchronized
% to the level you want.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatlabMPI
% Dr. Jeremy Kepner
% MIT Lincoln Laboratory
% kepner@ll.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Get size and rank.
  comm_size = MPI_Comm_size(comm);
  my_rank = MPI_Comm_rank(comm);

  if(comm_size < 2)
    disp('Too few processors (need at least 2) to synchronize');
    return;
  end


  % Create a unique tag id for the synch message.
  synch_tag = 99999;

  if (my_rank == starter_rank)
    % Get a zero clock.
    zero_clock = clock;

    % Get current time
    current_time = etime(clock,zero_clock);

    % Compute synchronized start time.
    start_time = current_time + delay;

    % Broadcast to everyone else.
    MPI_Bcast( starter_rank, synch_tag, comm, zero_clock, start_time );

    % Get current time again.
    current_time = etime(clock,zero_clock);

    % Compute pause time.
    pause_time = start_time - current_time

    % Pause.
%    pause(pause_time);

  else

    % Receive time from starter.
    [zero_clock start_time] = MPI_Recv( starter_rank, synch_tag, comm );

    % Get current time again.
    current_time = etime(clock,zero_clock);

    % Compute pause time.
    pause_time = start_time - current_time

    % Pause.
%    pause(pause_time);

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2002 Massachusetts Institute of Technology
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

