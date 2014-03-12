%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script implements a basic image convolution
% across multiple processors.  It illustrates how to
% use overlaping boundaries.
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pBlurimage
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pBlurimage
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pBlurimage',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set image size (scaled by numlabs), filter size and blur.
Nx = 2^11*Np;  Ny = 1024; Nk = 2.^5;  Nblur = 2;
%Nx = 2^9*Np;  % Debug.

PARALLEL = 1;  % Set control flag.
CHECK = 0;  % Check answer with serial calculation.

Zmap = 1;    % Create serial map.
if (PARALLEL)
  Zmap = map([Np 1],{},0:Np-1,[Nk 0]);  % Create parallel map with overlap.
end

% Create starting image and working images..
if (CHECK)   im = zeros(Nx,Ny); end


Z = zeros(Nx,Ny,Zmap) + 1.e-4;   % Create 2D distributed array.
[ii jj] = find(sprand(Nx,Ny,1e-4));  % Create non-zeros.
for i=1:numel(ii)
  Z(ii(i),jj(i))=1;    % Insert non-zeros.
end

myI = global_ind(Z,1);  % Get local i indices.
myJ = global_ind(Z,2);  % Get local j indices.

Z0 = Z;

if (CHECK)
  im = zeros(Nx,Ny) + 1.e-4;
  for i=1:numel(ii)
   im(ii(i),jj(i))=1;
  end
end

kernel = ones(Nk,Nk);  % Create kernel.

tic;
Z = synch(Z);  % Synchronize boundary conditions.
Tlaunch = toc;
disp(['Launch Time (sec)                  = ',num2str(Tlaunch)]);

tic;  % Set start time.

for iblur = 1:Nblur    % Loop over each blur.
  Zloc = local(Z);  % Get local data.
  Zloc(1:end-Nk+1,1:end-Nk+1) =conv2(Zloc,kernel,'valid');    % Perform covolution.
  Z = put_local(Z,Zloc);  % Put local back.
  Z = synch(Z);    % Copy overlaping boundaries.
  if (CHECK)
    im(1:end-Nk+1,1:end-Nk+1) = conv2(im,kernel,'valid');
  end

end

Tcompute = toc;  % Get blur time.
disp(['Compute Time (sec)                 = ',num2str(Tcompute)]);

% Compare results.
if (CHECK)
  maxDiff = max(max(abs(local(Z) - im(myI,myJ))));
  disp(['Max error                        = ',num2str(maxDiff)]);
end

totalOps = 2.*Nblur*Nk*Nk*Nx*Ny; % Number of operations.

% Print compute performance.
gigaflops = totalOps / Tcompute / 1.e9;
disp(['Performance (Gigaflops)            = ',num2str(gigaflops)]);

% Display on leader.
if (Pid == 0)
  figure
  plot(ii,jj,'o'); axis([0 Nx 0 Ny],'square');
  figure
  imagesc(rot90(local(Z),1));  axis('square','off');
%  save workspace
end
