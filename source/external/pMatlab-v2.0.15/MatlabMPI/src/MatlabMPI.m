% Message Passing Interface (MPI) for Matlab
%
%                 ------- Core Lite Profile -------         
%    MPI_Run              Runs a matlab script in parallel.
%    MPI_Init             Inititializes at the beginning.
%    MPI_Comm_size        Gets number of processors in a communicator.
%    MPI_Comm_rank        Gets rank of current processor within a communicator.
%    MPI_Send             Sends a message to a processor (non-blocking).
%    MPI_Recv             Receives message from a processor (blocking).
%    MPI_Finalize         Cleans up at the end.
%
%
%                 ------- Core Profile -------
%    MPI_Abort            Function to kill all matlab jobs
%                         started by MatlabMPI.
%    MPI_Bcast            Broadcast a message (blocking).
%    MPI_Probe            Returns a list of all incoming messages.
%    MPI_cc               Compiles using Matlab mcc.
%
%
%                 ------- Core Plus Profile -------
%
%                           [No functions, yet.]
%
%                 ------- User Utility functions -------
%
%    MatMPI_Host_rank
%
%    MatMPI_Delete_all    MatlabMPI function to delete all files created
%                         by MatlabMPI.
%    MatMPI_Comm_settings Can be edited by users to change the
%                         behavior of MatlabMPI (unix/windows, rsh/ssh, ...).
%    MatMPI_Save_messages MatlabMPI function to prevent messages
%                         from being deleted; useful for debugging purposes.
%    MatMPI_Comm_dir      MatlabMPI function for switching directory
%                         used for carrying out communication.
%
%                 ------- Library Utility functions -------
%
%    MatMPI_Buffer_file   MatlabMPI function for generating
%                         buffer file name.  Used by MPI_Send and MPI_Recv.
%    MatMPI_Lock_file     MatlabMPI function for generating
%                         lock file name.  Used by MPI_Send and MPI_Recv.
%    MatMPI_Commands      MatlabMPI function for generating
%                         commands.  Used by MPI_Run.
%    MatMPI_Comm_init     MatlabMPI function for creating MPI_COMM_WORLD.
%                         Used by MPI_Run.
%    MatMPI_Sleep
%
%    MatMPI_dir_map
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

