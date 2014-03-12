%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZoomImage: zoom in on an image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pZoomImage
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pZoomImage
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pZoomImage',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set image size, number frames, start and stop scale.
N = 512;  Ns  = 16; Sstart = 32; Send = 1;
N = 256; % Debug.
sigma  = 0.5;    % Width of blur kernel.
PARALLEL = 1;    % Set control flag.
Zmap = 1;       % Create serial map.
if (PARALLEL)  
  Zmap = map([1 1 Np],{},0:Np-1);  % Create parallel map.
%  dist(1).dist = 'b';  dist(2).dist = 'b'; dist(3).dist = 'c';
%  Zmap = map([1 1 Np],dist,0:Np-1);  % 1D Cyclic Distribution.
end
Z = zeros(N,N,Ns,Zmap);    % Create distributed array.
S = linspace(Sstart,Send,Ns);  % Compute scale factors.
disp('Zooming frames...');
tic;    % Start clock.
Z0 = referenceFrame(N,0.1,0.8);  % Create reference frame.
% Compute local frames.
Zloc = zoomFrames(Z0,S( global_ind(Z,3) ),sigma);
Tcompute = toc;   % Stop Clock.
disp(['Compute Time (sec)                 = ',num2str(Tcompute)]);
Z  = put_local(Z,Zloc);  % Insert into distributed  array.
tic;  % Start Clock.
Zagg = agg(Z);    % Aggregate on leader.
Tcomm = toc;   % Stop Clock.
disp(['Launch+Comm Time (sec)             = ',num2str(Tcomm)]);

% Compute gigaflops.
N_k = ceil(S*(5*sigma));
totalOps = 2.*sum((N_k.^2)).*(N.^2);
GigaFlops = 1.e-9*totalOps/Tcompute;
disp(['Performance (Gigaflops)            = ',num2str(GigaFlops)]);


% Display on leader.
if (Pid == 0)
  figure(1); clf;
  set(gcf,'Name','Simulated Frames','DoubleBuffer','on','NumberTitle','off');
  for frameIndex=[1:Ns]
    imagesc(squeeze(Zagg(:,:,frameIndex)));
    drawnow;
  end
end

%%%%%%%%%%%%%%%%%%%%%
