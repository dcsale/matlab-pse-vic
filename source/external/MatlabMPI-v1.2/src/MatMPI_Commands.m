function [defscommands, unix_command] = ...
  MatMPI_Commands(m_file,rank,MPI_COMM_WORLD)
% MatMPI_Commands  -  Commands to launch a matlab script remotely.
%
%  [defscommands, unix_command] = ...
%    MatMPI_Commands(m_file,rank,MPI_COMM_WORLD)
%

  % Unix vs. Windows file seperator.
            dir_sep = '/';
  if (ispc) dir_sep = '\'; end
            file_ext = '.sh';
  if (ispc) file_ext = '.bat'; end

  % Get host.
%  if (isunix) host = getenv('HOST');          end
  if (isunix) [status, host] = unix('hostname'); host = deblank(host); end
  if (ispc)   host = getenv('computername');  end

  qq = '"';
  sp = ' ';

  % Check if this is a compiled script.
  COMPILED_FLAG = 0;
  % Syntax is 'compiled m_file.exe'.
  compile_sep = findstr(m_file,'compiled ');  % Look for 'compiled'.
  exe_sep = findstr(m_file,'.exe');  % Look for '.exe'.
  if (compile_sep)
    if (exe_sep)
      COMPILED_FLAG = 1;
      % Parse m_file to recover m_file name.
      m_file = m_file(1,10:(exe_sep-1));
    end
  end

  % Set newline string.
  nl = sprintf('\n');

  % Create filename each Matlab job will run at startup.
  defsbase = ['MatMPI/MatMPIdefs' num2str(rank)];
  defsfile = [defsbase '.m'];
  comm_mat_file = 'MatMPI/MPI_COMM_WORLD.mat';
  outfile = ['MatMPI/' m_file '.' num2str(rank) '.out'];

  % Get single quote character.
  q = strrep(' '' ',' ','');

  % Create Matlab MPI setup commands.
  commands{1} = ['global MPI_COMM_WORLD;' nl];
  commands{2} = ['load ' q comm_mat_file q ';' nl];
  commands{3} = ['MPI_COMM_WORLD.rank = ' num2str(rank) ';' nl];
%  commands{4} = ['delete(' q defsfile q ');' nl];
  commands{5} = [m_file ';' nl];

  defscommands = '';

  % Get info on the target machine.
  machine_id = MPI_COMM_WORLD.machine_id(1,rank+1);
  machine = MPI_COMM_WORLD.machine_db.machine{1,machine_id};
  remote_launch = MPI_COMM_WORLD.machine_db.remote_launch{1,machine_id};
  remote_flags = MPI_COMM_WORLD.machine_db.remote_flags{1,machine_id};
  matlab_command = MPI_COMM_WORLD.machine_db.matlab_command{1,machine_id};
  type = MPI_COMM_WORLD.machine_db.type{1,machine_id};

  % Print name of the target machine we are launching on.
  disp(['Launching MPI rank: ' num2str(rank) ' on: ' machine]);

  % Create base matlab command.
  matlab_command = [matlab_command ' < ' defsfile ' > ' outfile ];
  if (strcmp(type,'pc'))  % Target is a pc.
    matlab_command = [ ...
      'matlab /nodesktop /minimize /nosplash /r MatMPIdefs' num2str(rank) ...
      ' /logfile ' outfile];
  end
  if (COMPILED_FLAG)
     % Run compiled script with rank as a command line argument.
     matlab_command = ['./' m_file '.exe ' num2str(rank) ' > ' outfile ];
  end
  % Here are some other versions I have tried that don't work well on unix.
  % matlab_command = [matlab_command ' -r ' defsfile ' -logfile ' outfile ' > /def/null'];
  % matlab_command = [matlab_command ' -r ' defsfile ' -logfile ' outfile ];
  % matlab_command = [matlab_command ' -r ' defsfile ' > ' outfile ];

  % Determine how to run script and where to send output.
  if (strcmp(machine,host))  % Target is host.
    % if (rank == 0)  
    % Check if running with host& set.
    if ((rank == 0) && (MatMPI_Host_rank(MPI_COMM_WORLD) == 0))
      % Run defsfile scipt interactively.
      defscommands = [commands{1} commands{2} commands{3} commands{5}];
      unix_command = nl;
    else
      if(not(COMPILED_FLAG))
        % Write commands to a .m text file.
        fid = fopen(defsfile,'wt');
        n_command = size(commands,2);
        for i_command=1:n_command
          fwrite(fid,commands{i_command});
        end
        fclose(fid);
      end
      % Create command to run defsfile locally and pipe output to another file.
      unix_command = [matlab_command ' &' nl 'touch MatMPI/pid.' machine '.$!' nl];
      if(strcmp(type,'pc'))  % Target is a pc.
%        unix_command = ['start /b ' matlab_command nl];
%        unix_command = ['start /b /high ' matlab_command nl];
%        unix_command = ['start /b ' matlab_command nl 'copy nul ' machine '.pc' nl];
        unix_command = ['start /b ' matlab_command nl 'copy nul MatMPI\pid.' machine '.pc' nl];
        % PC equivalent to touch is 'copy nul filename.txt'

      end
    end
  else  % Target is a remote machine.

    if(not(COMPILED_FLAG))
      % Write commands to a .m text file.
      fid = fopen(defsfile,'wt');
      n_command = size(commands,2);
      for i_command=1:n_command
        fwrite(fid,commands{i_command});
      end
      fclose(fid);
    end

    % Create command to run defsfile locally and pipe output to another file.
    unix_command = [matlab_command ' &' nl 'touch MatMPI/pid.' machine '.$!' nl];

    % Remote machine is a pc.
    if(strcmp(type,'pc'))  % Target is a pc.
%       unix_command = ['start /b ' matlab_command nl];
%       unix_command = ['start /b /high ' matlab_command nl];
%       unix_command = ['start /b ' matlab_command nl 'copy nul ' machine '.pc' nl];
       unix_command = ['start /b ' matlab_command nl 'copy nul MatMPI\pid.' machine '.pc' nl];
       % PC equivalent to touch is 'copy nul filename.txt'

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

