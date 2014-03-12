%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pBeamformer is an example of a simple parallel beamformer.
% The first part generates synthetic
% data and then passes the data into to the second half
% which processes the data.
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pBeamformer
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pBeamformer
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pBeamformer',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set number of time snapshots, sensors, beams and frequencies.
Nt = 600; Ns = 90;  Nb = 40; Nf = 100;
Nt = 100; % Debug.
PARALLEL = 1;  % Set control flag.
Xmap = 1;      % Create serial map.
if (PARALLEL)  
  Xmap = map([1 1 Np],{},0:Np-1);  % Create parallel map.
end

% ALLOCATE PARALLEL DATA STRUCTURES ---------------------

X0 = zeros(Nt,Nb,Nf,Xmap);   % Source array.
X1 = sqrt(Ns).*complex(rand(Nt,Ns,Nf,Xmap),rand(Nt,Ns,Nf,Xmap));  % Sensor input.
X2 = zeros(Nt,Nb,Nf,Xmap);  % Beamformed output.
X3 = zeros(Nt,Nb,Np,Xmap);  % Intermediate sum.
myI_f = global_ind(X1,3); % Get local indices.
% Get local parts of arrays.
X0loc = local(X0);  X1loc = local(X1);  X2loc = local(X2);

% CREATE STEERING VECTORS ---------------------

% Pick an arbitrary set of frequencies.
freq0 = 10;  frequencies = freq0 + (0:Nf-1);
% Create local steering vector by passing local frequencies.
myV = squeeze(Beamformer_vectors(Ns,Nb,frequencies(myI_f)));

% STEP 0: Insert targets ---------------------

% Insert two targets at different angles.
X0loc(:,round(0.25*Nb),:) = 1;  X0loc(:,round(0.5*Nb),:) = 1;


% STEP 1: CREATE SYNTHETIC DATA. ---------------------
tic;	% Start timer.

for i_t=1:Nt  % Loop over time snapshots.
  for i_f=1:numel(myI_f)   % Loop over local frequencies.
    % Convert from beams to sensors.
    X1loc(i_t,:,i_f) = X1loc(i_t,:,i_f).' + ...
      (squeeze(myV(:,:,i_f)) * squeeze(X0loc(i_t,:,i_f)).');
  end
end

% STEP 2: BEAMFORM AND SAVE DATA. ---------------------

for i_t=1:Nt   % Loop over time snapshots.
  for i_f=1:numel(myI_f)      % Loop over local frequencies.
    % Convert from sensors back to beams.
    X2loc(i_t,:,i_f) = ...
      abs(squeeze(X1loc(i_t,:,i_f)) * squeeze(myV(:,:,i_f))).^2;
  end
end

for i_f=1:numel(myI_f)   % Loop over frequencies.
   X_i_f = squeeze(X2loc(:,:,i_f));  % Get a matrix of data.
   filename = ['dat/pBeamformer_freq.' num2str(myI_f(i_f)) '.mat'];
   save(filename,'X_i_f');      % Save to a file.
end

Tcompute = toc;
disp(['Compute Time (sec)                 = ',num2str(Tcompute)]);

% STEP 3: SUM ACROSS FREQUNCY. ---------------------

X3 = put_local(X3,sum(X2loc,3));   % Sum local part and put into X3.
tic;  % Start timer.
x3 = squeeze(sum(agg(X3),3)); % Aggregate X3 on leader and complte sum.
Tcomm = toc;
disp(['Launch+Comm Time (sec)             = ',num2str(Tcomm)]);

% STEP 4: Finalize and display. ---------------------
% Display on leader.
if (Pid == 0)
  imagesc( abs(squeeze(X0loc(:,:,1))) );
  pause(1.0);
  figure;
  imagesc( abs(squeeze(X1loc(:,:,1))) );
  pause(1.0);
  figure;
  imagesc( abs(squeeze(X2loc(:,:,1))) );
  pause(1.0);
  figure;
  imagesc(x3);
end

%whos X0 X1 X2 X3 X0loc X1loc X2loc x3 myI_f myV
