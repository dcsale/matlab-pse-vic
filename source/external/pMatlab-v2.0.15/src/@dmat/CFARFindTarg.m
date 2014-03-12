function R = CFARFindTarg( C, T )
% check that the input is distributed along the correct axis
% recall C(Nbeam, Nrg, Ndop) should _not_ be broken up along
% Nrg.  Ideally, handle both cases where it is broken up across
% beams or dopplers, but for testing purposes remap distribution
% is anything other than across beams:
  
if C.dim == 2 && T.dim == 2
  % Check C
  g = grid(C);          % distribution grid scheme
  proc_list = g(:)';    % list of processors uses
  grid_dims = size(g);  % grid size
  old_map   = C.map;    % map method used
  dist_spec = old_map.dist_spec;  % distribution specification
  
  if grid_dims(2) ~= 1  % the matrix is broken up along the columns
    warning(['@dmat/CFAR: the C matrix is incorrectly mapped ' ...
             ', remapping along rows.']);
    grid_spec = [grid_dims(1)*grid_dims(2) 1];
    new_map = map(grid_spec, dist_spec, proc_list);
    C = remap(C, new_map);
  end
  
  % Check T
  g = grid(T);          % distribution grid scheme
  proc_list = g(:)';    % list of processors uses
  grid_dims = size(g);  % grid size
  old_map   = T.map;    % map method used
  dist_spec = old_map.dist_spec;  % distribution specification
  
  if grid_dims(2) ~= 1  % the matrix is broken up along the columns
    warning(['@dmat/CFAR: the T matrix is incorrectly mapped ' ...
             ', remapping along rows.']);
    grid_spec = [grid_dims(1)*grid_dims(2) 1];
    new_map = map(grid_spec, dist_spec, proc_list);
    T = remap(T, new_map);
  end  
  
  % Create the output map
  R = zeros( C.size(1), C.size(2), C.map);
  
  % Now operate on the local data
  local_size = localdims(C.falls, C.dim);
  local_C = local(C);
  local_T = local(T);
  local_R = local(R);
  
  local_R = (local_C.^2)./local_T;
  
  R = put_local(R, local_R);
  
else
  error('Incorrect matrix dimensions');
end
