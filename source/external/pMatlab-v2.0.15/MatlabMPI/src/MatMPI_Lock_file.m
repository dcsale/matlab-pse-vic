function lock_file = MatMPI_Lock_file(source,dest,tag,comm)
% MatMPI_Lock_file  -  function for creating lock file name.
%
%   lock_file = MatMPI_Lock_file(source,dest,tag,comm)
%

  machine_id = comm.machine_id(1,dest+1);
  dir = comm.machine_db.dir{1,machine_id};
             sep = '/';
  if (ispc)  sep = '\'; end
  lock_file = [dir,sep,'p',num2str(source),'_p',num2str(dest),'_t',num2str(tag),'_lock.mat'];
  
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


