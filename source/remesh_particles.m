function [xp, wp, PART] = remesh_particles(SIM, MESH, wf)
%% ideas
% could also use this function to creat particles inbetween a min and max
% value, for example, could be used to create particle isosurfaces !

%% =======================================================================%
% create the particle locations at the nodes within vorticity support, plus
% some additional padding (ghost particles?)
% ========================================================================%
wf_mag     = sqrt(wf{1}.^2 + wf{2}.^2 + wf{3}.^2);
tol        = 1e-4;                          % NOTE: should be function input and part of SIM structure
ii         = wf_mag > tol;                  % indicies where field is greater than tolerance
PART.nPart = sum(sum(sum(ii)));             % number of particles
PART.hp    = SIM.h_cutoff * MESH.dx(1);     % a smoothing radius (i.e., a cutoff length or core size)

xp = zeros(3, PART.nPart); % init the particles
type = 'collocated';
switch type
    case 'collocated'
        % create particles at the nodes
        xp(1,:) = MESH.xf(ii);
        xp(2,:) = MESH.yf(ii);
        xp(3,:) = MESH.zf(ii);
        wp      = interp_M2P(MESH, xp, wf);       % init weights by interpolating the field
    case 'staggered'
        % create particles at cell centers
%         wf_cen  = interp_M2M(SIM, MESH, MESH);
        xp(1,:) = MESH.xf_cen(ii);
        xp(2,:) = MESH.yf_cen(ii);
        xp(3,:) = MESH.zf_cen(ii);      
        wp      = interp_M2P(MESH, xp, wf_cen);   % init weights by interpolating the field
    otherwise
        error('unrecognized type of mesh')
end




