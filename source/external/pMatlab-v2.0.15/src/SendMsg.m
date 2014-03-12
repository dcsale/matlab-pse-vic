function SendMsg(dest, tag, varargin )
% Send message from current Matlab instance to dest.
  %globals
  global pMATLAB;

  if (numel(dest) == 1)
    MPI_Send(dest,tag,pMATLAB.comm,varargin);
  elseif (numel(dest) > 1)
    MPI_Mcast(Pid,dest,tag,pMATLAB.comm,varargin);
  end
