
function v = Beamformer_vectors(Nsensors,Nbeams,myFreqs)
%
% steering_vectors - routine to return broadband 3D focused single-path replica vector
%
% Input:       params.az = azimuth angles (deg), 0 = fwde       (1xM)
%              params.el = D/E angle (deg); 0= horizontal bea        (1xL)
%              focus_range = focus range (m) (OPTIONAL)                         (1x1)
%                defaults to a very large number to get plane-wave replica
%              params.arrayGeom = structure containing array element (x,y,z) location
%              params.freqs = vector of frequencies (Hz)                (1xF)
%              params.soundSpeed = speed of sound (m/sec)
% 
% Output:     v = steering vector (Nelements x M x L x F)
%


% Hard some code parameters.
params.el = 0;
params.az = linspace(0,360,Nbeams);

% frequencies and array positions are dimensionless.
params.freqs=myFreqs;
params.arrayGeom.x = linspace(-1000,1000,Nsensors);
params.arrayGeom.y = zeros(size(params.arrayGeom.x));
params.arrayGeom.z = zeros(size(params.arrayGeom.x));
params.numEls = length(params.arrayGeom.x);

params.soundSpeed=1500;


% get dimensions 
numElev = length(params.el);
numAz = length(params.az);
numFreqs = length(params.freqs);

% set so all azimuths have the same focus range
if isfield(params,'focus_range')==0,  % default to far-field
	focus_range = 1e10*ones(1,numAz);
else
	focus_range = params.focus_range*ones(1,numAz);
end


% shoehorn rr structure into P_array.
P_array=[params.arrayGeom.x' params.arrayGeom.y' params.arrayGeom.z'];
P_array_matrix=P_array(:,:,ones(1,numAz));

v = zeros(params.numEls,numAz,numElev,numFreqs);

for ielev = 1:numElev,

	% Define the vector that points at this azimuth and elevation 
	% from the array phase center
	pointing_vectors(1,:) = cos(params.az*pi/180)*cos(params.el(ielev)*pi/180);
	pointing_vectors(2,:) = sin(params.az*pi/180)*cos(params.el(ielev)*pi/180);
	pointing_vectors(3,:) = ones(1,numAz)*sin(params.el(ielev)*pi/180);

	% Compute the actual focus point (meters) 
	focus_points = pointing_vectors * diag(focus_range);

	% Compute the difference in range to each element in the array with 
	% respect to the array phase center
	focus_points_matrix = reshape(kron(focus_points,ones(params.numEls,1)),...
							params.numEls,3,numAz);

	delta_range = sqrt(squeeze(sum((P_array_matrix - ...
					focus_points_matrix).^2,2))) - ones(params.numEls,1)*focus_range;

	% Compute the true array response vectors to the source 
	% (azimuth,elevation,focus range)

	for ifrq = 1:numFreqs,
		freq = params.freqs(ifrq);
		v(:,:,ielev,ifrq) = exp(j*2*pi*delta_range*freq/params.soundSpeed);
	end


end
