function varargout = MPI_Recv( source, tag, comm )
% MPI_Recv  -  Receives message from source.
%
%  [var1, var2, ...] = MPI_Recv( source, tag, comm )
%
%    Receives message from source with a given tag
%    and returns the variables in the message.
%
%    source can be an iteger from 0 to comm_size-1
%    tag can be any integer
%    comm is an MPI Communicator (typically a copy of MPI_COMM_WORLD)
%

  % Get processor rank.
  my_rank = MPI_Comm_rank(comm);

  % Get file names.
  buffer_file = MatMPI_Buffer_file(source,my_rank,tag,comm);
  lock_file = MatMPI_Lock_file(source,my_rank,tag,comm);

  % Spin on lock file until it is created.
  loop = 0;
% exist() creates a lot of NFS traffice, better to use fopen.
%  while exist(lock_file) ~= 2
%    loop = loop + 1;
%  end

  temp_fid = fopen(lock_file,'r');

%i = 1;
  while temp_fid == -1
%disp(num2str(i));
    temp_fid = fopen(lock_file,'r');
%i = i + 1;
    loop = loop + 1;
    % Sleep statement allows cleaner profiling, but adds latency.
%    MatMPI_Sleep;
%    if (ispc)
%      MatMPI_Recv_sleep(0.1);
%    end
  end

  fclose(temp_fid);

  % Read all data out of buffer_file.
  buf = load(buffer_file);
  % Delete buffer and lock files.
  if (not(comm.save_message_flag))
    delete(buffer_file);
%    MatMPI_Sleep;
%    pause(0.1);
    delete(lock_file);
  end
  % Get variable out of buf.
  varargout = buf.varargin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatlabMPI
% Dr. Jeremy Kepner
% MIT Lincoln Laboratory
% kepner@ll.mit.edu
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

