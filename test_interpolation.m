

%% Preamble
close all;
fclose all;
profile off;
diary off;
clc;

% Add paths to the source code
if ~isdeployed
    % add folder to source code
    dir_src = pwd;
    addpath([dir_src filesep 'source'])
    addpath([dir_src filesep 'source' filesep 'PoissonFFT'])
    
    % add external dependencies
    dir_ext = [pwd filesep 'source' filesep 'external']; 
    addpath([dir_ext filesep 'bipolar_colormap'])
    addpath([dir_ext filesep 'consolidator'])
    addpath([dir_ext filesep 'interparc'])
    addpath([dir_ext filesep 'plot3k'])
    addpath([dir_ext filesep 'ODE_Solvers'])
    addpath([dir_ext filesep 'write_VTK' filesep 'vtk_writers'])

end

%% Test Accuracy of P2M and M2P interpolation

% read the input file
inputFile = 'ctrl.m';
[SIM, ENV, MESH, CTRL] = VortexParticle_Init(inputFile);    

% init the mesh
MESH = init_mesh(SIM, MESH);

% init a field on the mesh
wf = init_field(CTRL, SIM, MESH, ENV);

% initialize particle representation of field
[xp, wp, PART] = remesh_particles(SIM, MESH, wf);

% initialize particles at every mesh node, and then interpolate the field onto the particles
xp_tmp      = zeros(SIM.dim, numel(MESH.xf{1}(:))); 
xp_tmp(1,:) = MESH.xf{1}(:); % init the particles at mesh nodes including in the ghost layer
xp_tmp(2,:) = MESH.xf{2}(:);
xp_tmp(3,:) = MESH.xf{3}(:);
wp_tmp      = interp_M2P(SIM, MESH, xp_tmp, wf);

% only keep particles above some cutoff tolerance
wp_mag_tmp = sqrt( wp_tmp(1,:).^2 + wp_tmp(2,:).^2 + wp_tmp(3,:).^2);
ii         = wp_mag_tmp > SIM.tol_remesh;
nPart_2    = sum(ii);
xp_2       = zeros(SIM.dim, nPart_2); % init the particles within tolerance cutoff
wp_2       = zeros(SIM.dim, nPart_2); % init the particles within tolerance cutoff
xp_2(1,:)  = xp_tmp(1,ii);
xp_2(2,:)  = xp_tmp(2,ii);
xp_2(3,:)  = xp_tmp(3,ii);
wp_2(1,:)  = wp_tmp(1,ii);
wp_2(2,:)  = wp_tmp(2,ii);
wp_2(3,:)  = wp_tmp(3,ii);

% now we should get the original field back if we interpolate P2M (next measure the error on this)
wf_2 = interp_P2M(SIM, MESH, xp_2, wp_2);

%% Error Analysis

% compute the errors
err_wf_x = max(max(max( wf{1} - wf_2{1} )));
err_wf_y = max(max(max( wf{2} - wf_2{2} )));
err_wf_z = max(max(max( wf{3} - wf_2{3} )));

% visualize the fields
plot_field(  wf, MESH, 'Vorticity 1')
plot_field(wf_2, MESH, 'Vorticity 2')

% visualize the particles
plot_particles(  xp,   wp, MESH, 'Vorticity 1')
plot_particles(xp_2, wp_2, MESH, 'Vorticity 2')




