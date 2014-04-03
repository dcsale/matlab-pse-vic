function [dir_pc, dir_unix] = ...
  MatMPI_dir_map(machine_db,dir)
% MatMPI_dir_map  -  Takes care of pc/unix translation of pwd.
%
%  [dir_pc, dir_unix] = ...
%    MatMPI_dir(machine_db,dir)
%

  % Default is to do nothing.
  dir_pc = dir;
  dir_unix = dir;

  % Check if a directory mapping has been defined.
  % If so, convert directory names.
  if( isfield(machine_db,'pc_unix_dir_map') )
    pc_base = machine_db.pc_unix_dir_map{1};
    pc_n = length(pc_base);
    unix_base = machine_db.pc_unix_dir_map{2};
    unix_n = length(unix_base);

    % Check if pc_dir has a unix base.
    % remove from string.
    if (strncmp(dir_pc,unix_base,unix_n))
      % Swap bases.
      dir_pc = [pc_base dir_pc(unix_n+1:end)];
    end

    % Replace '/' with '\'.
    dir_pc = strrep(dir_pc,'/','\');

    % Check if unix_dir has a pc base.
    % remove from string.
    if (strncmpi(dir_unix,pc_base,pc_n))
      % Swap bases.
      dir_unix = [unix_base dir_unix(pc_n+1:end)];
    end

    % Replace '\' with '/'.
    dir_unix = strrep(dir_unix,'\','/');

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

