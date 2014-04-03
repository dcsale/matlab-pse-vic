function MPI_Abort()
% MPI_Abort  -  Aborts any currently running MatlabMPI sessions.
%
%   MPI_Abort()
%
%   Will abort any currently running MatlabMPI sessions.
%   by looking for leftover Matlab jobs and killing them.
%   Cannot be used after MatMPI_Delete_all.
%

  % Unix vs. Windows host name.
%  if (isunix) host = getenv('HOST');          end
  if (isunix) [status, host] = unix('hostname'); host = deblank(host); end
  if (ispc)   host = getenv('computername');  end

  % Get possibly user defined settings.
  machine_db_settings = MatMPI_Comm_settings;
  remote_launch = machine_db_settings.remote_launch;
  remote_flags = machine_db_settings.remote_flags;

  % Get list of pid files.
  pid_files = dir('MatMPI/pid.*.*');
  s = size(pid_files);
  n_files = s(1);

  % Create single qoute.
  q = strrep(' '' ',' ','');
  qq = '"';

  % Check if there are any files
  if (n_files < 1)
    disp('No pid files found');
  else

    % Loop over each file.
    for i_file=1:n_files

      % Get file name.
      file_name = deblank(pid_files(i_file).name);
      % Check if there is a pid appended.
      dir_sep = findstr(file_name,'.');
      if (dir_sep)
        
        % Parse file name.
        machine = file_name(1,(dir_sep(1)+1):(dir_sep(end)-1));
        pid = file_name(1,(dir_sep(end)+1):end);

        % Check if the target machine is a PC
        if (strcmp(pid, 'pc'))

           % Get the username
           if (isunix)
              [status, username] = unix('whoami');
              username = deblank(username);
           end
           if (ispc)
              username = getenv('USERNAME');
           end

           % Check if the target machine is the host
 	   if (strcmp(machine, host))
              % Don't do anything, otherwise MPI_Abort will kill the instance of Matlab that launched
              % the MPI_Abort command
              % unix_command = [ 'taskkill /f /FI "USERNAME eq ' username '" /im matlab.exe' q ];
              unix_command = '';
	   else
              % WinXP
              unix_command = [ remote_launch machine remote_flags ' ' q ...
                              'taskkill /f /FI "USERNAME eq ' username '" /im matlab.exe' q ];

              % Win2K
              % unix_command = [ remote_launch machine remote_flags ' ' q ...
              % 'kill -f matlab.exe' q];

              if (ispc)
                 unix_command = [ remote_launch machine remote_flags ' ' qq ...
                                 'taskkill /f /FI \"USERNAME eq ' username '\" /im matlab.exe' qq ];
              end
           end

        else

           % Check if the target machine is the host
  	   if (strcmp(machine, host))
              unix_command = [ 'kill -9 ' pid];
           else
              unix_command = [ remote_launch machine remote_flags ' ' ...
                               q 'kill -9 ' pid q];
              if (ispc)
                unix_command = [ remote_launch machine remote_flags ' ' ...
                                 qq 'kill -9 ' pid qq];
              end
           end
        end

        disp(unix_command);
        unix(unix_command);

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


