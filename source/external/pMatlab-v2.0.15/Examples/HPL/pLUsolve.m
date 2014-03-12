
function x = pLUsolve(A,b)
% Solves A x = b on a column distributed NxN matrix A
% using parallel LU factorization.

  [N,N] = size(A);                 % Get size of A.
  myJ = global_block_range(A,2);   % Get the local columns.

%tic;
  [L,U,pivots] = pLUfactor(A);     % Call parallel LU.
%TpLUfactor = toc
%tic;
  Lloc = local(L);                 % Get local L.
  Uloc = local(U);                 % Get local U.
  x = b(pivots,:);                 % Pivot rows of b.

  for p = 0:Np-1                   % Loop up over each Pid.
    tag = mod(p,32);               % Set message tag.
    if Pid == p        
      i = (myJ(2)+1):N;            % Upper block of rows.
      j = 1:(myJ(2)-myJ(1)+1);     % Lower block of columns.
      k = myJ(1):myJ(2);           % Middle block of rows.
      x(k,:) = Lloc(k,j)\x(k,:);   % Solve L x' = x.
      x(i,:) = x(i,:) - Lloc(i,j)*x(k,:);
      if p < Np-1
        SendMsg(p+1,tag,x)         % Send x to the higher Pid.
      end
    elseif Pid == p+1
      x = RecvMsg(p,tag);          % Recv x from the lower Pid.
    end
  end

  for p = (Np-1):-1:0              % Loop down over each Pid.
    tag = mod(p,32);               % Set message tag.
    if Pid == p                   
      i = 1:(myJ(1)-1);            % Upper block of rows.
      j = 1:(myJ(2)-myJ(1)+1);     % Lower block of columns.
      k = myJ(1):myJ(2);           % Middle block of rows.
      x(k,:) = Uloc(k,j)\x(k,:);   % Solve L x' = x.
      x(i,:) = x(i,:) - Uloc(i,j)*x(k,:);
      if p > 0
        SendMsg(p-1,tag,x);        % Send x to the lower Pid.
      end
    elseif Pid == p-1
      x = RecvMsg(p,tag);          % Recv x from the higher Pid.
    end
  end

  tag = 3;                         % Reset message tag.
  if Pid == 0;                    
    p_higher = 1:Np-1;
    if not(isempty(p_higher))
%      SendMsg(p_higher,tag,x);     % Send x to all other Pids.
      for pp = p_higher
        SendMsg(pp,tag,x); % Send Lloc and pPivots to higher Pids.
      end
    end
  else
    x = RecvMsg(0,tag);            % Received x from lower Pids.
  end
%TpLUsolve = toc

end
