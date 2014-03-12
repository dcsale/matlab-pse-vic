function MPI_COMM_WORLD = MatMPI_Comm_init(n_proc,machines)
% MatMPI_Comm_init  -  Creates generic communicator.
%
%   MPI_COMM_WORLD = MatMPI_Comm_init(n_proc,machines)
%


  % Get number of machines to launch on.
  n_machines = size(machines,2);
  n_m = max(n_machines,1);

  % Set default machine.
%  if (isunix) host = getenv('HOST');          end
  if (isunix) [status, host] = unix('hostname'); host = deblank(host); end
  if (ispc)   host = getenv('computername');  end
  machine = host;


  % Initialize comm.
  MPI_COMM_WORLD.rank = -1;
  MPI_COMM_WORLD.size = n_proc;
  MPI_COMM_WORLD.save_message_flag = 0;
  MPI_COMM_WORLD.group = (1:n_proc)-1;
  MPI_COMM_WORLD.machine_id = zeros(1,n_proc);
  MPI_COMM_WORLD.host_rank = 0;

  % Initialize machine database.
  machine_db.n_machine = n_m;			% Number of machines.
  machine_db.type = cell(1,n_m);		% 'unix' or 'pc'.
  machine_db.machine = cell(1,n_m);		% Machine names.
  machine_db.dir = cell(1,n_m);			% Communication directory.
  machine_db.matlab_command = cell(1,n_m);	% Matlab command.
  machine_db.remote_launch = cell(1,n_m);	% Remote launch command.
  machine_db.remote_flags = cell(1,n_m);	% Remote launch flags.
  machine_db.n_proc = zeros(1,n_m);		% # processes on this machine.
  machine_db.id_start = zeros(1,n_m);		% Start index.
  machine_db.id_stop = zeros(1,n_m);		% Stop index.

  % Start setting up machine id.
  for i_rank=0:n_proc-1
    i_machine = mod(i_rank,n_m) + 1;
    machine_db.n_proc(1,i_machine) = machine_db.n_proc(1,i_machine) + 1;
  end

  % Get possibly user settings.
  machine_db_settings = MatMPI_Comm_settings;
  if (isfield(machine_db_settings,'pc_unix_dir_map'))
    machine_db.pc_unix_dir_map = machine_db_settings.pc_unix_dir_map;
  end

  % Set default type
  default_type = 'unix';
  if (ispc) default_type = 'pc'; end

  % Set paths.
  [pwd_pc pwd_unix] = MatMPI_dir_map(machine_db,pwd);

  % Set machine_db values.
  for ii=1:n_m
    machine_db.type{1,ii} = default_type;
    machine_db.machine{1,ii} = host;
    machine_db.dir{1,ii} = [pwd '/MatMPI'];
    machine_db.matlab_command{1,ii} = machine_db_settings.matlab_command;
    machine_db.remote_launch{1,ii} = machine_db_settings.remote_launch;
    machine_db.remote_flags{1,ii} = machine_db_settings.remote_flags;
    if (ii == 1) 
      machine_db.id_start(1,ii) = 1;
      machine_db.id_stop(1,ii) = machine_db.id_start(1,ii) + machine_db.n_proc(1,ii) -1;
    else
      machine_db.id_start(1,ii) = machine_db.id_stop(1,ii-1) + 1;
      machine_db.id_stop(1,ii) = machine_db.id_start(1,ii) + machine_db.n_proc(1,ii) -1;
    end

    id_start = machine_db.id_start(1,ii);
    id_stop = machine_db.id_stop(1,ii);

    MPI_COMM_WORLD.machine_id(1,id_start:id_stop) = ii;

    % Check if there is a machines list.
    if (n_machines > 0)
      machine = machines{ii};
      machine_db.machine{1,ii} = machine;

      % Strip out '&' if present.
      amp_sep = findstr(machine,'&');
      % remove from string.
      if (amp_sep);
         machine = [machine(1:amp_sep-1) machine(amp_sep+1:end)];
         machine_db.machine{1,ii} = machine;
      end

      % Check if same as host.  DOESN'T HANDLE host:dir syntax.
      if (strcmp(machine,host))
        % Set type to type of host.
        machine_db.type{1,ii} = default_type;
      else
        % Use user specified default (probably 'unix').
        machine_db.type{1,ii} = machine_db_settings.type;
      end

      % Check if ':unix' or ':pc' is appended, if so, override type.
      unix_sep = findstr(machine,':unix');
      if (unix_sep);
        machine_db.type{1,ii} = 'unix';
        % remove from string.
        machine = [machine(1:unix_sep-1) machine(unix_sep+5:end)];
        machine_db.machine{1,ii} = machine;
      end
      pc_sep = findstr(machine,':pc');
      if (pc_sep)
        machine_db.type{1,ii} = 'pc';
        % remove from string.
        machine = [machine(1:pc_sep-1) machine(pc_sep+3:end)];
        machine_db.machine{1,ii} = machine;
      end

      % Change default directory according to type.
      type = machine_db.type{1,ii};
      if (strcmp(type,'pc'))
        machine_db.dir{1,ii} = [pwd_pc '\MatMPI'];
      end
      if (strcmp(type,'unix'))
        machine_db.dir{1,ii} = [pwd_unix '/MatMPI'];
      end

      % Check if there is a directory appended.
      dir_sep = findstr(machine,':');
      if (dir_sep)
        dir_piece = machine(1,(dir_sep(1)+1):end);
        machine = machine(1,1:dir_sep(1)-1);
        machine_db.dir{1,ii} = dir_piece;
        machine_db.machine{1,ii} = machine;
      end

      if (strcmp(machine,host))
        if (amp_sep)
          MPI_COMM_WORLD.host_rank = -1;
        end
      end

    end
  end

  % Add machine_db to communicator.
  MPI_COMM_WORLD.machine_db = machine_db;

  % Write out.
  comm_mat_file = 'MatMPI/MPI_COMM_WORLD.mat';
  save(comm_mat_file,'MPI_COMM_WORLD');

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

