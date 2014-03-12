
function [L,U,pivots] = pLUfactor(A)
% Parallel LU factorization on a column distributed NxN matrix A.
  global pMATLAB;

  [N,N] = size(A);                          % Get size of A.  
  myJ  = global_block_range(A,2);           % Get the local columns.
  allJ = global_block_ranges(A,2);          % Get all the local columns.
  allNloc = allJ(:,3) - allJ(:,2) + 1;      % Number of local columns.
  Aloc = local(A);                          % Get local part of A.
  Nloc = size(Aloc,2);                      % Get local column size.
  A = put_local(A,0);                       % Zero out local part of A.

  % Turn off message deletion to allow non-blocking broadcasts.
  pMATLAB.comm = MatMPI_Save_messages(pMATLAB.comm,1);

  pivots = (1:N)';                          % Initialize pivots.

  for p = 0:Np-1                            % Loop over each Pid.
    disp(num2str(p));
    tagHigh = mod(2*p+1,32);                % High message tag.
    tagLow = mod(2*p,32);                   % Low message tag.

    if p == Pid                             % Factor block p.


      i = myJ(1):N;                         % Get rows to work on.

if 0
%  tic;
      AlocSub = Aloc(myJ(1):N,:);
      [AlocSub pPivots] = dgetrf(AlocSub);
      Lloc = tril(AlocSub);
      Lloc(eye(numel(i),Nloc)==1) = 1;      % Set diagonal to 1.
%  Tgetrf = toc
end


if 1
%  tic;

      AlocSub = Aloc(myJ(1):N,:);
      [AlocSub Uloc pPivots] = lu(AlocSub,'vector');  % Factor.

      AlocSub(1:Nloc,:) = AlocSub(1:Nloc,:) + (Uloc - eye(Nloc,Nloc));    % Add back upper part.

      Lloc = tril(AlocSub);               % Get lower triangle.
      Lloc(eye(numel(i),Nloc)==1) = 1;      % Set diagonal to 1.

%  Tlu = toc
end
      pPivots = pPivots+myJ(1)-1;           % Update pivots.


%tic;
      pHigh = (p+1):Np-1;                   % Higher Pids.
      if not(isempty(pHigh))                
%       SendMsg(pHigh,tagHigh,Lloc,pPivots); % Send Lloc and pPivots to higher Pids.
        for pp = pHigh
          SendMsg(pp,tagHigh,Lloc,pPivots); % Send Lloc and pPivots to higher Pids.
        end
      end

      pLow = 0:(p-1);                       % Lower Pids.
      if not(isempty(pLow))                 
%        SendMsg(pLow,tagLow, pPivots);      % Send pPivots to lower Pids.
        for pp = pLow
          SendMsg(pp,tagLow, pPivots);      % Send pPivots to lower Pids.
        end
      end
%Tsend = toc
      Aloc(myJ(1):N,:) = AlocSub;

    elseif Pid > p
%tic;
      [Lloc pPivots] = RecvMsg(p,tagHigh);  % Receive L and pivots.
%Trecv = toc
    elseif Pid < p
%tic;
      pPivots = RecvMsg(p,tagLow);          % Receive pivots.
%Trecv = toc
    end

%tic;
    i = allJ(p+1,2):N;                      % Remaining rows.
    pivots(i) = pivots(pPivots);            % Update pivot list.

    if p ~= Pid
      Aloc(i,:) = Aloc(pPivots,:);          % Apply pivots everywhere else.
    end

    % Apply multipliers to following blocks.
    if Pid > p
      iL1 = 1:allNloc(p+1);                 % Upper rows of Lloc.
      iL2 = (allNloc(p+1)+1):size(Lloc,1);  % Lower rows of Lloc.
      iA1 = allJ(p+1,2):allJ(p+1,3);        % Current p rows od Aloc.
      iA2 = (allJ(p+1,3)+1):N;              % Lower rows of Aloc.
      Aloc(iA1,:) = Lloc(iL1,:)\Aloc(iA1,:);   % Solve.
      Aloc(iA2,:) = Aloc(iA2,:) - Lloc(iL2,:)*Aloc(iA1,:);  % Update.
    end

%Tapply = toc

  end

%tic;
  % Turn  message deletion back on.
  % This shouldn't be necessary, but it is a good habit.
  pMATLAB.comm = MatMPI_Save_messages(pMATLAB.comm,0);

  Lloc = tril( Aloc, -(myJ(1)) );           % Get lower triangle.
  i = myJ(1):myJ(2);                        % Get local rows.
  Lloc(i,:) = Lloc(i,:) + speye(Nloc,Nloc); % Add ones to diagonal.
  L = put_local(A,Lloc);                    % Put into distributed array.

  Uloc = triu( Aloc, -(myJ(1)-1) );         % Get upper triangle.
  U = put_local(A,Uloc);                    % Insert into U.

%Tput = toc

end
