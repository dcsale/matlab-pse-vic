function MPI_Send( dest, tag, comm, varargin )
% MPI_Send  -  Sends variables to dest.
%
%  MPI_Send( dest, tag, comm, var1, var2, ...)
%
%    Send message containing variables to dest with a given tag
%
%    dest can be an iteger from 0 to comm_size-1
%    tag can be any integer
%    comm is an MPI Communicator (typically a copy of MPI_COMM_WORLD)
%

  % Get processor rank.
  my_rank = MPI_Comm_rank(comm);

  % Create buffer and lock file.
  buffer_file = MatMPI_Buffer_file(my_rank,dest,tag,comm);
  lock_file = MatMPI_Lock_file(my_rank,dest,tag,comm);

  % Save buf to file.
  % Save using Matlab 6 (R13) MAT file format
%  if (str2num(version('-release')) > 13)
   if exist('OCTAVE_VERSION','builtin')
     save(buffer_file,'varargin');
   else
     save(buffer_file,'varargin','-v6');
   end
%  else
%    save(buffer_file,'varargin');
%  end

  % Create lock file.
%  fclose(fopen(lock_file,'w'));
  fclose(fopen(lock_file,'w+'));
%  [fid,errmsg] = fopen(lock_file,'w+');
%  fid
%  errmsg
%  if (fid == -1)
%    errmsg
%  end
%  fclose(fid);


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

