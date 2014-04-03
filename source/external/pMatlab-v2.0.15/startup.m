%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example startup.m file user can copy into
%   ~/matlab/startup.m (Linux & Mac)
%   C:\My Documents\MATLAB (Windows)
%
%   These suggestions may not work with your environment.
%   Please make sure that
%   Matlab can pick up the startup.m when it starts.
%
%   Please start Matlab from the command line in 
%   order to run multiple Matlab processes on a
%   single computer.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add MatlabMPI.
addpath('/home/kepner/pMatlab/MatlabMPI/src');

% Add pMatlab.
addpath('/home/kepner/pMatlab/src');

% To launch with PCT uncomment the following line.
%addpath('/home/kepner/pMatlab/PCTstubs');
%Type 'help pRUN' to see how to use.

if exist('OCTAVE_VERSION','builtin')
  addpath('/home/kepner/matlab');

  % Octave users should also call their
  % startup.m file in their .octaverc

  % Use --no-window-system in MPI_Run
  % or set (0, 'defaultfigurevisible', 'off');

end

% Windows users uncomment to add .\MatMPI to path:
% if ispc
%addpath .\MatMPI
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER CONFIGURATION PARAMETERS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rehash;  % Rehash everything.

% Note: if you type "clear all" at the command line you
% will need to rerun this function.
pMatlabGlobalsInit;

% Put everything below here in pMatlabGlobalsInit?

if not(exist('OCTAVE_VERSION','builtin'))
  % R2009b says maxNumCompThreads will be deprecated in the future release
  % Comment this function out.
  % maxNumCompThreads(1);
end
