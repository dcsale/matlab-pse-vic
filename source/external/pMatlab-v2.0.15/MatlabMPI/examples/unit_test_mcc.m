function status = unit_test_mcc( m_file, n_proc, machines, time_limit, time_check )
%
% Test a matlab program.
%

  % Launch the script.
  disp(['Running: ' m_file]);
  MPI_Run(MPI_cc(m_file), n_proc, machines);

  % Determine how many times to do checks.
  n_check = round(time_limit / time_check);

  % Check status.
  i_check = 1;
  keep_checking = 1;
  for i=1:n_check
    pause(time_check);
    % Flush buffers.
    unix(' ');
    % Check to see how many are done.
    [s,n_d] = unix(['grep SUCCESS ./MatMPI/' m_file '.*.out | wc -l']);
    n_done = str2num(n_d);
    disp([num2str(time_check*i) ' seconds.  Completed: ' num2str(n_done)]);
    if (n_done == n_proc)
      status = 'SUCCESS';
      disp([m_file ': SUCCESS']);
      MatMPI_Delete_all;
      pause(10.0);
      return;
    end
  end    

  status = 'FAIL';
  disp([m_file ': FAIL']);
  MPI_Abort;
  unix(['cp -r MatMPI ' m_file '.' num2str(n_proc) '.FAIL.MatMPI']);
  MatMPI_Delete_all;
  pause(10.0);

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

