function T = CFARNoiseEst( C, Ncfar, G )
% Pattern Match Algorithm implemented in pMatlab as part of development
% process for pMapper implementation.  See pMapper/main/pCFAR.m for 
% detailed description of the algorithm
  
% check that the input is distributed along the correct axis
% recall C(Nbeam, Nrg, Ndop) should _not_ be broken up along
% Nrg.  Ideally, handle both cases where it is broken up across
% beams or dopplers, but for testing purposes remap distribution
% is anything other than across beams:
  
if C.dim == 2
  g = grid(C);          % distribution grid scheme
  proc_list = g(:)';    % list of processors uses
  grid_dims = size(g);  % grid size
  old_map   = C.map;    % map method used
  dist_spec = old_map.dist_spec;  % distribution specification
  
  % Check that C is broken up along beams; if not, remap
  if grid_dims(2) ~= 1  % the matrix is broken up along the columns
    warning(['@dmat/CFAR: the matrix is not mapped along the appropriate ' ...
             'dimension, remapping along rows.']);
    grid_spec = [grid_dims(1)*grid_dims(2) 1];
    new_map = map(grid_spec, dist_spec, proc_list);
    C = remap(C, new_map);
  end
  
  % Create the output map
  T = zeros( C.size(1), C.size(2), C.map);
  
  % Now operate on the local data
  local_size = localdims(C.falls, C.dim);
  local_C = local(C);
  local_T = local(T);
  
  % We only work on squared variables in C
  local_C = local_C .^ 2;
  
  for ii=1:local_size(1)
    for jj=1:local_size(2)
      % If this is the first range gate, compute using brute force.
      % Assumes there are enough range gates to compute off the 
      % right hand size of the sliding window.  MATLAB will generate
      % an index out of bounds error if this is false
      if jj == 1;         
        local_T(ii,jj) = sum( local_C(ii,jj+G+1:jj+Ncfar) ) / (2 * Ncfar);
      else
        if jj+1+G+Ncfar <= local_size(2)
          local_T(ii,jj) = local_T(ii,jj) + local_C(ii, jj+1+G+Ncfar); 
        end
        if jj-G >= 1
          local_T(ii,jj) = local_T(ii,jj) + local_C(ii, jj-G);
        end
        if jj-G-Ncfar >= 1
          local_T(ii,jj) = local_T(ii,jj) - local_C(ii, jj-G-Ncfar);
        end
        if jj+G+1 <= local_size(2)
          local_T(ii,jj) = local_T(ii,jj) - local_C(ii, jj+G+1);
        end
        local_T(ii,jj) = (local_T(ii,jj) / 2 * Ncfar) + local_T(ii,jj-1);
      end
    end
  end
  
  T = put_local(T, local_T);
  
else
  error('Incorrect matrix dimensions');
end
