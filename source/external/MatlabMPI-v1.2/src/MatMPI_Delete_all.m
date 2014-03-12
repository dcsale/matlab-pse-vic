function MatMPI_Delete_all()
% MatMPI_Delete_all  -  Deletes leftover MatlabMPI files.
% 
%  MatMPI_Delete_all()
%
%

  dir_sep = '/';
  if (ispc) dir_sep = '\'; end

  % Check MatMPI directory exists.                                              
  if exist('./MatMPI', 'dir') == 0                                                
    disp('Nothing to delete.');                                                 
    return                                                                      
  end                                                                           

  % Check MPI_COMM_WORLD exists.                                              
  if exist('./MatMPI/MPI_COMM_WORLD.mat') == 0
    disp('No MPI_COMM_WORLD, deleting anyway, files may be leftover.');
  else

    % First load MPI_COMM_WORLD.
    load './MatMPI/MPI_COMM_WORLD.mat';

    % Set newline string.
    nl = sprintf('\n');
    % Get single quote character.
    q = strrep(' '' ',' ','');

    % Get number of machines.
    n_m = MPI_COMM_WORLD.machine_db.n_machine;

    % Loop backwards over each machine.
    for i_m=n_m:-1:1

      % Get number of processes to launch on this machine.
      n_proc_i_m = MPI_COMM_WORLD.machine_db.n_proc(1,i_m);

      if (n_proc_i_m >= 1)

        % Get communication directory.
        comm_dir = MPI_COMM_WORLD.machine_db.dir{1,i_m};

        % Delete buffer and lock files in this directory.
        delete([comm_dir dir_sep 'p*_p*_t*_buffer.mat']);
        delete([comm_dir dir_sep 'p*_p*_t*_lock.mat']);
      end
    end
  end

  % Delete MatMPI directory.
 % rmdir(['.' dir_sep 'MatMPI' dir_sep '*']);
  rmdir(['.' dir_sep 'MatMPI'],'s');

%  if (isunix)
%    delete(['.' dir_sep 'MatMPI']);
%  else
%    rmdir(['.' dir_sep 'MatMPI']);
%    delete(['.' dir_sep 'MatMPI']);
%    dos(['rmdir /S /Q .' dir_sep 'MatMPI']);
%  end

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

