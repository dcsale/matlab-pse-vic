function [xp, wp, PART] = remesh_particles(SIM, MESH, wf)
%% ideas
% could also use this function to creat particles inbetween a min and max
% value, for example, could be used to create particle isosurfaces

%% update the mesh to have tighter bounds on the vorticity support (giving better particle resolution)
% MESH = updateMesh(SIM, MESH, xp_old, PART);

%% =======================================================================%
% create the particle locations at the nodes within vorticity support
% NOTE: should we also add some additional padding (ghost particles?)
% ========================================================================%
wf_mag     = sqrt(wf{1}.^2 + wf{2}.^2 + wf{3}.^2);
ii         = wf_mag > SIM.tol_remesh;      	% indicies where field is greater than tolerance


PART.nPart = sum(sum(sum(ii)));            	% number of particles
PART.hp    = SIM.h_cutoff * max(MESH.dx); 	% a smoothing radius (i.e., a cutoff length or core size)
PART.volp  = prod(MESH.dx);                 % particle volumes (using midpoint/rectangular rule)

%% some error checking 
if PART.nPart == 0 
    error('[Error] no particles created during remesh.')
end

% crucial to check particle overlap condition hp/dx > 1
if (PART.hp / max(MESH.dx)) < 1.0
    error('[ERROR] particle overlap condition is violated! hp/dx > 1')
end

%%
xp = zeros(SIM.dim, PART.nPart); % init the particles
type = 'collocated';
switch type
    case 'collocated'
        % create particles at the nodes (this might include the ghost nodes)
        xp(1,:) = MESH.xf{1}(ii);
        xp(2,:) = MESH.xf{2}(ii);
        xp(3,:) = MESH.xf{3}(ii);
        wp      = interp_M2P(SIM, MESH, xp, wf);       % init weights by interpolating the field
    case 'staggered'
        % create particles at cell centers (not finished)
%         wf_cen  = interp_M2M(SIM, MESH, MESH);
%         xp(1,:) = MESH.xf_cen(ii);
%         xp(2,:) = MESH.yf_cen(ii);
%         xp(3,:) = MESH.zf_cen(ii);      
%         wp      = interp_M2P(SIM, MESH, xp, wf_cen);   % init weights by interpolating the field
    otherwise
        error('[ERROR] unrecognized type of mesh')
end

