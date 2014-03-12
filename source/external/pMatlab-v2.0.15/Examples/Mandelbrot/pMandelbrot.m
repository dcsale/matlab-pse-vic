%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Mandelbrot set in parallel.
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pMandelbrot
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pMandelbrot
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pMandelbrot',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set number of iterations, mesh size and threshold.
N=2400;  Niter=20;  epsilon = 0.001;
N=800;  % Debug
PARALLEL = 1;                   % Set control flag.
Wmap = 1;                       % Create serial map.
if (PARALLEL)
  Wmap = map([Np 1],{},0:Np-1);        % Create parallel map.
%  dist(1).dist = 'c';  dist(2).dist = 'b';
%  Wmap = map([Np 1],dist,0:Np-1);  % 1D Cyclic Distribution.
end
W = zeros(N,N,Wmap);            % Create distributed array
Wloc = local(W);                % Get local part.
myI = global_ind(W,1);  % Get local i indices.
myJ = global_ind(W,2);  % Get local j indices.
[ReC ImC] = meshgrid(myJ./(N/2) -1.6, myI./(N/2) -1 );
Cloc = complex(ReC,ImC);        % Initialize C.
Zloc = Cloc;                    % Initialize Z.
ieps = 1:numel(Wloc);           % Initialize indices.
tic;                            % Start clock.
for i=1:Niter;                  %  Compute Mandelbrot set.
  Zloc(ieps) = Zloc(ieps) .* Zloc(ieps)  + Cloc(ieps);
  Wloc(ieps) = exp(-abs(Zloc(ieps)));
  ieps = ieps( find(Wloc(ieps) > epsilon) );
end
W = put_local(W,Wloc);         % Put back into W;
Tcompute = toc;               % Stop clock.
disp(['Compute Time (sec)                 = ',num2str(Tcompute)]);

tic;                            % Start Clock.
W1 = agg(W);                    % Aggregate back to leader.
Tcomm = toc;         % Stop clock.
disp(['Launch+Comm Time (sec)             = ',num2str(Tcomm)]);

% Display on leader.
if (Pid == 0)
%  colormap copper(256);
%  pcolor(W1);
%  shading flat;
  imagesc(W1);
  axis('square','equal','off');
end

%%%%%%%%%%%%%%%%%%%%%
