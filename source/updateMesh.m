function MESH = updateMesh(MESH, pad, xp, PART)
%% =======================================================================%
% update the mesh: 
% ========================================================================%
% parent mesh
if MESH.adaptive
    MESH.x  = linspace(min(xp(1,:)) - pad, max(xp(1,:)) + pad, MESH.Nx)'; % NOTE: once variable core sizes are implemented PART.hp will be a particle quantity, so query the core size of the particle at the mesh boundaries
    MESH.y  = linspace(min(xp(2,:)) - pad, max(xp(2,:)) + pad, MESH.Ny)';
    MESH.z  = linspace(min(xp(3,:)) - pad, max(xp(3,:)) + pad, MESH.Nz)';
    MESH.dx = MESH.x(2) - MESH.x(1);
    MESH.dy = MESH.y(2) - MESH.y(1);
    MESH.dz = MESH.z(2) - MESH.z(1);
    [MESH.X, MESH.Y, MESH.Z] = meshgrid(MESH.x, MESH.y, MESH.z);
end

% crucial to check hp/dx > 1
if ~(PART.hp / PART.h_cutoff) > 1.0
    error('[ERROR] condition hp/dx > 1 is violated!')
end

end % updateMesh()

