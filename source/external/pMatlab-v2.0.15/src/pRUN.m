function defscommands = pRUN( m_file, n_proc, machines )
% pRUN  -  Run m_file on multiple processors.
% Supports the following modes.
%  PARALLEL = 0;
%    eval(pRUN('pZoomImage',1,{}));
%  PARALLEL = 1;
%    eval(pRUN('pZoomImage',1,{}));
%    eval(pRUN('pZoomImage',2,{}));
%    eval(pRUN('pZoomImage',2,{'cluster'}));
%    eval(pRUN('pZoomImage',4,{'cluster'}));
%
% If using Mathworks PCT, machines = {} will run locally,
% otherwise machines needs to be a matlab scheduler object
% typically obtained from the findResource command.

%  if (n_proc < 2)
%    % Just return the script.
%    defscommands = m_file;
%  elseif (n_proc > 1)
  if (n_proc > 0)
    disp(['Submitting ' m_file ' on ' num2str(n_proc) ' processor(s).']);

    % Set wrapper file name.
    pRUN_Parallel_Wrapper_file = 'pRUN_Parallel_Wrapper';

      script_file_head = '.pRUN_Parallel_Stub_';
      script_file_tail = '_temp';

      script_file = [script_file_head,m_file,script_file_tail];

      % Run locally in pmode.
      % Abort left over jobs.
      MPI_Abort;
      pause(2.0);

      % Delete left over MPI directory
      MatMPI_Delete_all;
      delete([script_file_head,'*',script_file_tail]);
      pause(2.0);

%      fclose(fopen(script_file,'w+t'));
      [fid,errmsg] = fopen(script_file,'w+t');
      fclose(fid);
      defscommands = MPI_Run(pRUN_Parallel_Wrapper_file, n_proc, machines);

    % Use has decided to use PCT to launch program.
    pctWrapper = 'pctRUN_Parallel_Wrapper';
    if not(isempty(which(pctWrapper)))

      global pctJOB;

      pctconfig('hostname', 'localhost');

      % Create a parallel job.
      sched = findResource('scheduler','type','local');
      % Check if machines is a scheduler object, if not, run locally.
      [temp,classname,temp] = fileparts(class(machines));
      if (strcmp(classname,'distcomp'))
        sched = machines;
      end

      pctJOB = createParallelJob(sched);
      set(pctJOB, 'MaximumNumberOfWorkers', n_proc-1);
      set(pctJOB, 'MinimumNumberOfWorkers', n_proc-1);

      % Create task to call eval on my_file;
      tsk = createTask(pctJOB, 'eval', 0, {pctWrapper});
      set(tsk, 'CaptureCommandWindowOutput', true);

      % Submit the job.
      submit(pctJOB);
      waitForState(pctJOB,'running');
      %pause(2.0);

      % Add pctFinalize to locally run commands.
      defscommands = [defscommands ' pctFinalize;'];

    end
  end

  % Prepend pMatlabGlobalsInit in case it was cleared.
  defscommands = ['pMatlabGlobalsInit; ' defscommands ];

return
end
