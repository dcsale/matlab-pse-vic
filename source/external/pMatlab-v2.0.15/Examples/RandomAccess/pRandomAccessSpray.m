%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core of RandomAccess benchmark using Spray communication pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ran = ranStarts;  % Initialize random sequence.
tag = 0; % Initialize message tag.

% Set optimal send and receive order to minimize communication contention.
mySendOrder = circshift([0:(Pid-1) (Pid+1):(Np-1)],[0, -Pid]);
myRecvOrder = fliplr(circshift([0:(Pid-1) (Pid+1):(Np-1)], [0, -Pid]));

% Loop over all update blocks.
for ib = myBLOCK
%  disp(['Update block #' num2str(ib)]);
  tag = mod(tag + 1,32);   % Increment messsage tag.

  % Create random numbers using official HPC Challenge RandomAccess
  % random number generator
  ran = RandomAccessRand(ran);
  Xi = double(bitand(ran,mask))+1;   % Compute global table index.

  for p = mySendOrder               % Find and send updates.
    ranSend = ran( (Xi >= allX(p+1,2)) & (Xi <= allX(p+1,3)) );
    SendMsg(p,tag,ranSend);
  end

  % Get local updates.
  ranRecv = ran( (Xi >= myX(1)) & (Xi <= myX(2)) );

  for p = myRecvOrder       % Receive updates.
    ranRecv = [ranRecv RecvMsg(p,tag)];  % Append receives.
  end

  Xi = double(bitand(ranRecv,mask))+1-(myX(1)-1);  % Compute local index.

  if (not(VALIDATE))     % Fast update.
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
