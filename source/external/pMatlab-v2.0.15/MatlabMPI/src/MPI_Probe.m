function [message_rank, numeric_tag, string_tag] = MPI_Probe(source, tag, comm)
% MPI_Probe  -  Returns a list of all messages waiting to be received.
% 
%  [message_rank, numeric_tag] = MPI_Probe(source, tag, comm)
%
%    Source can be an integer or a wildcard '*'.
%    Tag can be an integer or a string, including wildcards.
%    Message_rank is a column-vector of source ranks, which may be empty.
%    Numeric_tag is a corresponding column-vector of numeric tags;
%        messages with non-numeric tags are not reported.
%
%  [message_rank, numeric_tag, string_tag] = MPI_Probe(source, tag, comm)
%
%    Numeric_tag is a corresponding column-vector of numeric tags;
%        non-numeric tags are reported as NaN.
%    String_tag is a column-cellarray of strings, both numeric or
%        non-numeric.

  % Get processor rank.
  my_rank = MPI_Comm_rank(comm);

  % Get lock file names.
  lock_file = MatMPI_Lock_file(source, my_rank, tag, comm);

  % Get and count the pending messages.
  message_files = dir(lock_file);
  n_files = length(message_files);

  % Create arrays to store rank and tag.
  rank = cell(n_files,1);
  string_tag = rank;

  % Parse out the source rank and tag for each message.
  rt = regexp({message_files.name}, '^p(\d+)_p\d+_t(.*)_lock', 'tokens');

  % Extract the rand and tag character strings.
  for msg = 1:n_files
    [rank{msg} string_tag{msg}] = deal(rt{msg}{1}{:});
  end

  % Convert to numeric.
  message_rank = str2double(rank);
  numeric_tag = str2double(string_tag);

  % Omit non-numeric tags if string_tags are not requested.
  if nargout == 2
      numericTags = ~isnan(numeric_tag);
      numeric_tag = numeric_tag(numericTags);
      message_rank = message_rank(numericTags);
  end

  % Check if there are no files
  if isempty(message_rank)
    % Sleep statement allows cleaner profiling, but adds latency.
    MatMPI_Sleep;
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


