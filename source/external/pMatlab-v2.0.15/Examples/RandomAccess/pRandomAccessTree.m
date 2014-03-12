%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core of RandomAccess benchmark using Tree communication pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ran = ranStarts;  % Initialize random sequence.
tag = 0; % Initialize message tag.

tree = factor(Np); % Factor processors into a tree.

% Loop over all update blocks.
for ib = myBLOCK
%  disp(['Update block #' num2str(ib)]);
  tag = mod(tag + 1,32);   % Increment messsage tag.

  % Create random numbers using official HPC Challenge RandomAccess
  % random number generator
  ran = RandomAccessRand(ran);

  ranRecv = ran;  % Init input.

  % Initialize values for computing splits.
  midP = Np ./ tree(1);   modP = Np;   k = 0;

  % Loop over nodes of tree.
%  for i=0:(numel(tree)-1)
  while (modP > 1)

    % Compute Pid to exchange info with.
    pairPid = mod((Pid + midP),modP) + floor(Pid./modP).*modP;

    % Compute Pid for splitting data.
    splitPid = midP + floor(Pid./modP).*modP;

    % Compute indices of all values.
    Xi = double(bitand(ranRecv,mask)) + 1;

    % Find values above and below split.
    hi = ranRecv( Xi >= allX(splitPid+1,2) );
    lo = ranRecv( Xi <  allX(splitPid+1,2) );

    % Exchange hi/lo data.
    if (Pid < pairPid)
      SendMsg(pairPid,tag,hi);
      ranRecv = [lo RecvMsg(pairPid,tag)];
    elseif (Pid > pairPid)
      SendMsg(pairPid,tag,lo);
      ranRecv = [hi RecvMsg(pairPid,tag)];
    end

    % Update values for computing splits.
    midP = midP./tree(k+1);   modP = modP./tree(k+1);   k=k+1;

  end

  % Select low bits for local table index.
  % Convert to double/add one for Matlab/Fortran style indexing.
  Xi = double(bitand(ranRecv,mask)) + 1 - (myX(1)-1);

  if (not(VALIDATE)) % Fast update.
    Xloc(Xi) = bitxor(Xloc(Xi),ranRecv);
  else    % Slow error free update.
    for j=1:length(ranRecv) 
      Xloc(Xi(j)) = bitxor(Xloc(Xi(j)), ranRecv(j));
    end
  end

end

% Insert local values back into global table.
if (VALIDATE)
  X = put_local(X,Xloc);
end
