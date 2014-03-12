function x = patternmatch(d, pattern)
%PATTERN_MATCH Pattern match algorithm on a distributed matrix.
%   patternmatch(Data, Pattern) uses the mean squared error between
%   each row of D and pattern (must be of the same length) to score
%   each entry.  The lower the score, the better the match.
%   If the matrix d is column distributed, patternmatch throws a
%   warning and redistributes the array.
%
%   Note: returned x is a new distributed array.  If d is an m x n
%   matrix, x is an m x 1 matrix.
  
if d.dim == 2  %distributed array is a matrix
  g = grid(d);                   % Distribution Grid Scheme
  proc_list = g(:)';             % Processors used
  grid_dims = size(g);           % Grid Size
  old_map = d.map;               % Map method used
  dist_spec = old_map.dist_spec; % Distribution Spec.
  
  % Check that the matrix is broken up along the rows, 
  % if not, remap
  if grid_dims(2) ~= 1   % the matrix is broken up along the columns
    warning(['@dmat/patternmatch: The matrix is not mapped along the ' ...
             'appropriate dimension, remapping along rows.']);
    grid_spec = [grid_dims(1)*grid_dims(2) 1];      % change to rows
    new_map = map(grid_spec, dist_spec, proc_list); % create new map
    d = remap(d, new_map);                          % remap
  end
  
  % Create the output distributed array
  x = zeros( d.size(1), 1, d.map );
  
  % Iterate across the local rows
  local_size = localdims(d.falls, d.dim);
  for ii = 1:local_size(1) 
    x.local(ii) = sum( ( d.local(ii) - pattern ).^2 ) / length(pattern);
  end
else
  error(['@dmat/patternmatch: patternmatch can only be applied to ' ...
         'distributed, 2D data matricies and non-distributed patterns']);
end
  
  
