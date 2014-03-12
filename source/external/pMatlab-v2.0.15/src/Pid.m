function x = Pid()
% Returns the Pid of the current Matlab instance.
  global pMATLAB;
  x = pMATLAB.my_rank;
end 
