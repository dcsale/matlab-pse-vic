function varargout = RecvMsg( source, tag )
% Receive message from source onto current Matlab instance.
  %globals
  global pMATLAB;

  varargout = MPI_Recv( source, tag, pMATLAB.comm );
