function [xp, PART] = init_particle_pos(PART, MESH, wf)
%% =======================================================================%
% create the particle locations at the nodes within vorticity support, plus
% some additional padding (ghost particles?)
% ========================================================================%
wf_mag     = sqrt(wf{1}.^2 + wf{2}.^2 + wf{3}.^2);
tol        = 1e-4; % NOTE: should be function input or part of SIM structure
ii         = wf_mag > tol;
PART.nPart = sum(sum(sum(ii)));

xp         = zeros(3, PART.nPart);

% create partiles at the nodes
xp(1,:)    = MESH.xf(ii);
xp(2,:)    = MESH.yf(ii);
xp(3,:)    = MESH.zf(ii);
% create partiles at the cell centers
% xp(1,:)    = MESH.xf(ii);
% xp(2,:)    = MESH.yf(ii);
% xp(3,:)    = MESH.zf(ii);

end % init_particle_pos
