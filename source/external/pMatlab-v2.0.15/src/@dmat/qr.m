function RQ = qr(A)
%% 
% RQ = qr(A) returns the QR decomposition of matrix A
%
%    Based on Algorithms 5.2.2 (Givens QR), 5.1.5, (Givens
%    Rotation) and 5.1.6 (Row Rotation) of:
%    "Matrix Computations" (2nd Ed) - Golub & Van Loan
%
%    Creates an orthonormal matrix Q and an upper-right
%    triangular matrix R such that A = Q'*R.  Givens
%    rotations are applied to zero out targeted 
%
%    The same permutation matricies are applied to the 
%    identity matrix to compute Q.
%
%    If A is an n-by-n matrix, RQ is an n-by-(2*n) matrix with
%    the upper-right triangular matrix R in the left half and 
%    Q the orthonormal matrix in the right half:
%
%    R(:,:) = RQ(1:n,1:n);
%    Q(:,:) = RQ(1:n,n+1:2*n);
%    
%    The algorithm operates on cols of the input matrix, and is 
%    most efficient when the distributed input matrix is row-wise
%    block distributed (otherwise the algorithm redistributes the
%    input dmat before processing)
%
%    The algorithm starts at the bottom of each column.  If 
%    the current column under processing is 'k', the algorithm
%    will zero out all A(rows > k, k).  Due to the distribution,
%    the upper-most row of each local matrix cannot be zeroed 
%    locally.  Thus, after each processor's local processing 
%    completes, the upper-most row of each local matrix must be
%    aggregated to a single processor and zeroed out by 
%    rotating it against row(k).  When this "global processing"
%    is completed, the modified rows are redistributed to the 
%    appropriate processors and local processing resumes on the 
%    next column k+1.
%%
  
global pMATLAB;

nThisProc = pMATLAB.my_rank;
nLastProc = pMATLAB.pList(end);

% Check on dimensions/distribution first

% Make sure it is 2-dimensional
if A.dim == 2
  gA = grid(A);
  grid_dimsA = size(gA);
  
  % Make sure input is row-distributed
  if grid_dimsA(2) ~= 1
    warning(['@dmat/qr: The input matrix is not mapped along' ...
             'the appropriate dimension, remapping along rows.']);
    
    grid_spec = [grid_dimsA(1)*grid_dimsA(2) 1];
    old_map = A.map;
    dist_spec = old_map.dist_spec;
    proc_list = gA(:)';
    new_map = map(grid_spec, dist_spec, proc_list);
    A = remap(A, new_map);
    gA = grid(A);
    grid_dimsA = size(gA);
  end

  % Create Outputs
  [Am An] = size(A);         % Global Size
  
  RQ = zeros( A.size(1), 2*A.size(2), A.map );
  
  % If it is square and rows are evenly dividable among the 
  % number of processors used for distribution:

  if (Am == An) && (rem(Am, length(pMATLAB.pList)) == 0)
 
    % Assume  no overlap
    [lAm lAn] = size(A.local); % Local Size (note, An == lAn)

    % Append appropriate portion of the identity matrix to 
    % the local portion of A  (Recall rank is zero-indexed)
    % lLZ = "left zeros" of identity matrix
    % lDO = "diagonal ones" of identity matrix
    % lRZ = "right zeros" of identity matrix
    lLZ = zeros(lAm, (nThisProc*lAm));
    lDO = diag(ones(1,lAm));
    lRZ = zeros(lAm, (nLastProc-nThisProc)*lAm);
    lAI = [A.local lLZ lDO lRZ];
    [lAIm lAIn] = size(lAI);

    
    % Begin processing
    for j = 1:An
      % Begin Local Processing
      %
      % Local processing has three possible code-paths to reflect
      % which part of the upper-right triangularization the
      % current column (j) is in for this processor (nThisProc).  
      % Recalling that the processor list is zero-indexed,
      % this could be one of:
      %
      %     "destined to become all zeros" ->
      %             j <= nThisProc * nLocalRows;
      %
      %     "destined to include the diagonal" ->
      %             j > nThisProc * nLocalRows &&
      %             j <= (nThisProc + 1) * nLocalRows
      %
      %     "to the right of the diagonal" ->
      %             j > (nThisProc + 1) * nLocalRows
      if j <= nThisProc*lAm;     % (destined to be all zeros)
        lower_bound = 2;         % zero to the top row, the
                                 % top row will be zeroed by
                                 % global processing
        
      elseif j <= (nThisProc+1)*lAm;  % includes the diagonal
        lower_bound = j - (nThisProc*lAm) + 1;
        
      else                       % to the right of the diag
        continue;
      end
      
      % Use Givens Rotations to zero out targeted indicies
      for i = lAm:-1:lower_bound
        %% Givens Algorithm, computes the rotation matrix 'P'
        %% (See beginning of this file for math references)
        if lAI(i,j) == 0    % already zero, don't rotate
          c = 1;
          s = 1;
        else
          if abs(lAI(i,j)) > abs(lAI(i-1,j))
            tau = -lAI(i-1,j) / lAI(i,j);
            s = 1 / sqrt( 1 + tau^2 );
            c = s * tau;
          else
            tau = -lAI(i,j) / lAI(i-1,j);
            c = 1 / sqrt( 1 + tau^2 );
            s = c * tau;
          end
        end
        
        P = [c s; -s c];
    
        %% Row Rotation, applies the rotation matrix R
        lAI(i-1:i,j:lAIn) = P' * lAI(i-1:i,j:lAIn);


      end
  
      % Global Processing - handled on the leader processor
      %
      % This should be a "barrier" - all processes come to 
      % this point, then wait for all other processes to arrive 
      % Then a single processor (the one that owns row j when
      % processing on column j, or the "pivot") will receive the top
      % (non-zeroed) row from processors that are attempting to 
      % zero out all rows in the active column (k) and, using
      % Givens rotations, zero each one out against the row 
      % that should be the non-zero A(k,k) on the diagonal.  
      % After all appropriate rows are zeroed, the updates are
      % communicated back to the appropriate processor, the
      % current row becomes k+1 and local processing resumes.
      
      % figure out if there is any global computation required
      % (none required for the last An/nProcs columns)
      nGlobalRotations = floor((An-j)/lAm);
      if( nGlobalRotations )
    
        % figure out who has the pivot
        nHasPivotProc = floor((j-1)/lAm);
    
        % only the processor with the pivot and processors 
        % after it in the processor list are active.  The processor
        % with the pivot receives, processes, and sends.  The remaining
        % active processors send, wait, receive.
        if( nThisProc > nHasPivotProc )
          nTag = (j * 10000) + (nThisProc * 10) + 0;
          MPI_Send( nHasPivotProc, nTag, pMATLAB.comm, lAI(1,j:lAIn) );
          lAI(1,j:lAIn) = MPI_Recv( nHasPivotProc, nTag+1, pMATLAB.comm );
        elseif( nThisProc == nHasPivotProc )
          for nSender = pMATLAB.pList(nThisProc+2:end)
            nTag = (j * 10000) + (nSender * 10) + 0;
            clear tMat;
            tMat(1,:) = lAI(j-(nThisProc*lAm),j:lAIn);
            tMat(2,:) = MPI_Recv( nSender, nTag, pMATLAB.comm );
            
            %% Givens Algorithm, computes the rotation matrix 'P'
            %% (See beginning of this file for math references)
            if tMat(2,1) == 0    % already zero, don't rotate
              c = 1;
              s = 1;
            else
              if abs(tMat(2,1)) > abs(tMat(1,1))
                tau = -tMat(1,1) / tMat(2,1);
                s = 1 / sqrt( 1 + tau^2 );
                c = s * tau;
              else
                tau = -tMat(2,1) / tMat(1,1);
                c = 1 / sqrt( 1 + tau^2 );
                s = c * tau;
              end
            end
        
            P = [c s; -s c];
    
            %% Row Rotation, applies the rotation matrix R
            tMat = P' * tMat;
            
            lAI(j-(nThisProc*lAm),j:lAIn) = tMat(1,:);
            MPI_Send( nSender, nTag+1, pMATLAB.comm, tMat(2,:) );
          end
        else
          ;; % do nothing, this processor is done
        end
      end
    end
    
    RQ = put_local(RQ, lAI);

  else %% if (square) and (size / nProcs is whole integer)
    error(['@dmat/qr: Input matrix must be square and ' ...
          'evenly dividable by the number of processors ' ...
          'to be distributed across']);
  end
else % not a 2D matrix
  error('@dmat/qr: Input must be a 2-dimensional matrix');
end
    
