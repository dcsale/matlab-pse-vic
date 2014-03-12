function runcommand = MPI_cc( m_file )
% MPI_cc  -  Compile a m_file for running with MPI_Run
%
%   runcommand = MPI_cc( m_file )
%
% To run, type:
%  
%    eval( MPI_Run( MPI_cc( m_file ), numprocesors, machines) );
%
% Or, type:
%
%    MPI_cc( m_file );
%    eval( MPI_Run( 'compiled m_file.exe' ), numprocesors, machines) );
%


  % Find wrappers directory.
  MPI_cc_dir = which('MPI_cc');
  index = findstr(MPI_cc_dir,'/MPI_cc.m');
  MPI_cc_dir = MPI_cc_dir(1:index-1);

  % Set directories and filenames.
  wrapper_dir = [MPI_cc_dir '/MatMPI_mcc_wrappers'];
  temp_dir = ['./.temp_MatMPI_mcc'];
  m_wrapper = 'my_function_wrapper.m';
  c_wrapper = 'my_function_p.c';
  defsfile = 'my_function.m';

  % Set newline string.
  nl = sprintf('\n');
  % Get single quote character.
  q = strrep(' '' ',' ','');

  % Create temporary directory.
  mkdir(temp_dir);

  % Prepend function header to user's script.
  unix(['cat ' wrapper_dir '/' m_wrapper ' ' m_file '.m  > ' temp_dir '/' defsfile]);

  % Create compilation command.
  mcc_command = ['mcc -v -t -W lib:multpkg -T link:exe -h ' temp_dir '/' defsfile ' ' wrapper_dir '/' c_wrapper ' libmmfile.mlib -d ' temp_dir ' -o ../' m_file '.exe']

  % Run command.
  eval(mcc_command);

  % Delete temp directory.
  rmdir(temp_dir,'s');

  % Return command that can be parsed by MPI_Run.
  runcommand = ['compiled ' m_file '.exe'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatlabMPI
% Dr. Jeremy Kepner
% MIT Lincoln Laboratory
% kepner@ll.mit.edu
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

