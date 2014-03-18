function MESH = updateMesh(SIM, MESH, xp, PART)
%% =======================================================================%
% update the mesh: account for periodic BC on particles, or growing the
% mesh if particle cross boundaries
% ========================================================================%
pad = SIM.pad * PART.hp;

% parent mesh
if MESH.adaptive
    % mesh boundaries
    MESH.xmin = [min(xp(1,:)) - pad, min(xp(2,:)) - pad, min(xp(3,:)) - pad]; 
    MESH.xmax = [max(xp(1,:)) + pad, max(xp(2,:)) + pad, max(xp(3,:)) + pad];
    % mesh coordinates
    MESH.x{1} = linspace(MESH.xmin(1), MESH.xmax(1), MESH.NX(1));
    MESH.y{2} = linspace(MESH.xmin(2), MESH.xmax(2), MESH.NX(2));
    MESH.z{3} = linspace(MESH.xmin(3), MESH.xmax(3), MESH.NX(3));
    % mesh spacing
    MESH.dx(1) = MESH.x{1}(2) - MESH.x{1}(1);
    MESH.dx(2) = MESH.x{2}(2) - MESH.x{2}(1);
    MESH.dx(3) = MESH.x{3}(2) - MESH.x{3}(1);
    % now add the ghost layer
    MESH.x{1} = linspace(MESH.xmin(1)-MESH.dx(1)*SIM.mbc, MESH.xmax(1)+MESH.dx(1)*SIM.mbc, MESH.NX(1)+2*SIM.mbc);
    MESH.x{2} = linspace(MESH.xmin(2)-MESH.dx(2)*SIM.mbc, MESH.xmax(2)+MESH.dx(2)*SIM.mbc, MESH.NX(2)+2*SIM.mbc);
    MESH.x{3} = linspace(MESH.xmin(3)-MESH.dx(3)*SIM.mbc, MESH.xmax(3)+MESH.dx(3)*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
    % mesh field
    [MESH.xf, MESH.yf, MESH.zf] = ndgrid(MESH.x{1}, MESH.x{2}, MESH.x{3});
    % centred domain
    centre         = (MESH.xmin + MESH.xmax)/2;
    MESH.xf_cen{1} = MESH.xf - centre(1);
    MESH.xf_cen{2} = MESH.yf - centre(2);
    MESH.xf_cen{3} = MESH.zf - centre(3); 
    
    % also need to update the particle volumes, but these evolve according
    % to the flow so must also solve additional system of ODE for particle
    % volumes  - IM NOT SURE ABOUT THIS
end

% crucial to check hp/dx > 1
if (PART.hp / max(MESH.dx)) < 1.0
    error('[ERROR] particle overlap condition is violated! hp/dx > 1!')
end

end % updateMesh()

