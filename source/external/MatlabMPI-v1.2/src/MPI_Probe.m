function [message_rank, message_tag] = MPI_Probe( source, tag, comm )
% MPI_Probe  -  Returns a list of all messages waiting to be received.
% 
%  [message_rank, message_tag] = MPI_Probe( source, tag, comm )
%
%    Source and tag can be an integer or a wildcard '*'.
%

  % Get processor rank.
  my_rank = MPI_Comm_rank(comm);

  % Get lock file names.
  lock_file = MatMPI_Lock_file(source,my_rank,tag,comm);

  % Check to see if there are any messages.
  message_files = dir(lock_file);
  n_files = length(message_files);

  % Create single qoute.
  q = strrep(' '' ',' ','');

  % Check if there are any files
  if (n_files < 1)
    % Set default (negative) return values.
    message_rank = '';
    message_tag = '';

    % Sleep statement allows cleaner profiling, but adds latency.
    MatMPI_Sleep;
  else
    % Create arrays to store rank and tag.
    message_rank = zeros(n_files,1);
    message_tag = message_rank;

    % Step through each file name and strip out rank and tag.
    for i_file=1:n_files

      % Get file name.
      file_name = message_files(i_file).name;

      % Parse the file name.
      [pathstr, name, ext, versn] = fileparts(file_name);
      if (length(name) > 0)
         [strings, count, errmsg, nextindex] = sscanf(name, 'p%d_p%d_t%d_lock');
      end

      % Sanity check; make sure that the file name was in the correct format
      if (size(strings, 1) > 0)
         message_rank(i_file) = strings(1);
         message_tag(i_file) = strings(3);
      end

    end

  end

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


