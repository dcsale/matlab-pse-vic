function x = Np()
% Returns the number Matlab instance currently running.
  global pMATLAB;
  x = pMATLAB.comm_size;
end
