global pctJOB;

waitForState(pctJOB,'finished');
outputmessages = get(pctJOB.Tasks, 'CommandWindowOutput');
errormessages = get(pctJOB.Tasks, 'ErrorMessage');

% Set wrapper file name.
pRUN_Parallel_Wrapper_file = 'pRUN_Parallel_Wrapper';

if (get(pctJOB,'MaximumNumberOfWorkers') == 1)
  dir_sep = '/';
  if (ispc) dir_sep = '\'; end

  output_fname = ['MatMPI' dir_sep  pRUN_Parallel_Wrapper_file ...
            '.' num2str(1) '.out'];
  fid = fopen(output_fname,'w+t');
  fwrite(fid,[outputmessages errormessages]);
  fclose(fid);
elseif (get(pctJOB,'MaximumNumberOfWorkers') > 1)
 for i_proc = 1:get(pctJOB,'MaximumNumberOfWorkers')
   output_fname = ['MatMPI/'  pRUN_Parallel_Wrapper_file ...
            '.' num2str(i_proc) '.out'];
   fid = fopen(output_fname,'w+t');
   fwrite(fid,[outputmessages{i_proc} errormessages{i_proc}]);
   fclose(fid);
  end
end

destroy(pctJOB);
clear('pctJOB','outputmessages','errormessages', ...
       'pRUN_Parallel_Wrapper_file','i_proc','fid','output_fname','dir_sep');
