% [fieldk, fieldx_cen, fieldx_ext, fieldr_ext, NX, vort] = init_field(NXs, Sim.testcase, plot_InitField);
function Mesh = init_mesh(CTRL, Sim, Mesh)

%--------------------------------------------------------------------------
%- Setup domain range and length
%--------------------------------------------------------------------------
if Sim.testcase == 1   
    Mesh.xmin = [-1, -1, -1]; 
    Mesh.xmax = [ 1,  1, 1];
elseif Sim.testcase == 2   
    Mesh.xmin = [0, 0, 0]; 
    Mesh.xmax = [2, 2, 2];
end
% L = xmax-xmin; % Domain length

%----------------------------------------------------------------------
% Setup domain (_ext is the extended domain)
%----------------------------------------------------------------------
% mesh coordinates
Mesh.x{1} = linspace(Mesh.xmin(1), Mesh.xmax(1), Mesh.NX(1));
Mesh.x{2} = linspace(Mesh.xmin(2), Mesh.xmax(2), Mesh.NX(2));
Mesh.x{3} = linspace(Mesh.xmin(3), Mesh.xmax(3), Mesh.NX(3));

% mesh spacing
Mesh.dx(1) = Mesh.x{1}(2) - Mesh.x{1}(1);
Mesh.dx(2) = Mesh.x{2}(2) - Mesh.x{2}(1);
Mesh.dx(3) = Mesh.x{3}(2) - Mesh.x{3}(1);

% mesh field
% [Mesh.x, Mesh.y, Mesh.z] = ndgrid(x{1}, x{2}, x{3});
[Mesh.xf, Mesh.yf, Mesh.zf] = ndgrid(Mesh.x{1}, Mesh.x{2}, Mesh.x{3});

% Centred domain
centre        = (Mesh.xmin + Mesh.xmax)/2;
Mesh.xf_cen{1} = Mesh.xf - centre(1);
Mesh.xf_cen{2} = Mesh.yf - centre(2);
Mesh.xf_cen{3} = Mesh.zf - centre(3);

% % Extended domain (FFT shifted)
% x_ext{1} = Mesh.x{1} - Mesh.xmin(1); 
% x_ext{1} = [x_ext{1} -x_ext{1}(end)-Mesh.dx(1) -x_ext{1}(end:-1:2)];
% 
% x_ext{2} = Mesh.x{2} - Mesh.xmin(2); 
% x_ext{2} = [x_ext{2} -x_ext{2}(end)-Mesh.dx(2) -x_ext{2}(end:-1:2)];
% 
% x_ext{3} = Mesh.x{3} - Mesh.xmin(3); 
% x_ext{3} = [x_ext{3} -x_ext{3}(end)-Mesh.dx(3) -x_ext{3}(end:-1:2)];
% 
% [Mesh.xf_ext{1}, Mesh.xf_ext{2}, Mesh.xf_ext{3}] = ndgrid(x_ext{1},x_ext{2},x_ext{3});
% Mesh.rf_ext = sqrt(Mesh.xf_ext{1}.^2 + Mesh.xf_ext{2}.^2 + Mesh.xf_ext{3}.^2);
% 
% if Sim.solve_vel == 2
%     % Wavenumbers for spectral differentiating
%     ks       = 1/Mesh.dx(1);
%     k_ext{1} = 2*pi*linspace(-ks/2,ks/2,2*Mesh.NX(1)+1);
%     k_ext{1} = fftshift(k_ext{1}(1:end-1));
% 
%     ks       = 1/Mesh.dx(2);
%     k_ext{2} = 2*pi*linspace(-ks/2,ks/2,2*Mesh.NX(2)+1);
%     k_ext{2} = fftshift(k_ext{2}(1:end-1));
% 
%     ks       = 1/Mesh.dx(3);
%     k_ext{3} = 2*pi*linspace(-ks/2,ks/2,2*Mesh.NX(3)+1);
%     k_ext{3} = fftshift(k_ext{3}(1:end-1));
% 
%     [Mesh.kf{1}, Mesh.kf{2}, Mesh.kf{3}] = ndgrid(k_ext{1}, k_ext{2}, k_ext{3});
% end

%     clear Mesh.x Mesh.y Mesh.z
        
end % function

% %% OLD WAY
% %% =======================================================================%
% % Set Initial Conditions
% % ========================================================================%
% switch SIM.example
%     case 'VortexRings'
%         % parent mesh
%         pad       = RING.Rmajor;
%         xMin      = min(RING.center_x) - pad*(max(RING.Rmajor) + max(RING.Rminor));
%         xMax      = max(RING.center_x) + pad*(max(RING.Rmajor) + max(RING.Rminor));
%         yMin      = min(RING.center_y) - pad*(max(RING.Rmajor) + max(RING.Rminor));
%         yMax      = max(RING.center_y) + pad*(max(RING.Rmajor) + max(RING.Rminor));
%         zMin      = min(RING.center_z) - pad*(max(RING.Rmajor) + max(RING.Rminor));
%         zMax      = max(RING.center_z) + pad*(max(RING.Rmajor) + max(RING.Rminor));
%            
%     case 'Turbine'
%         pad       = LL.ROTOR_DIA/2;
%         xMin      = -pad;
%         xMax      = +pad;
%         yMin      = -pad;
%         yMax      = +pad;
%         zMin      = 0;
%         zMax      = LL.HUB_HT + LL.ROTOR_DIA/2 + pad;
% 
%     otherwise
%         
% end
% %% parent mesh
% MESH.x        = linspace(xMin, xMax, MESH.NX(1))';
% MESH.y        = linspace(yMin, yMax, MESH.NX(2))';
% MESH.z        = linspace(zMin, zMax, MESH.NX(3))';
% MESH.dx       = MESH.x(2) - MESH.x(1);
% MESH.dy       = MESH.y(2) - MESH.y(1);
% MESH.dz       = MESH.z(2) - MESH.z(1);
% [MESH.X, MESH.Y, MESH.Z] = meshgrid(MESH.x, MESH.y, MESH.z);
% %% sub-grid scale mesh
% MESH_SGS.x    = linspace(xMin, xMax, MESH_SGS.NX(1))';
% MESH_SGS.y    = linspace(yMin, yMax, MESH_SGS.NX(2))';
% MESH_SGS.z    = linspace(zMin, zMax, MESH_SGS.NX(3))';
% MESH_SGS.dx   = MESH_SGS.x(2) - MESH_SGS.x(1);
% MESH_SGS.dy   = MESH_SGS.y(2) - MESH_SGS.y(1);
% MESH_SGS.dz   = MESH_SGS.z(2) - MESH_SGS.z(1);
% [MESH_SGS.X, MESH_SGS.Y, MESH_SGS.Z] = meshgrid(MESH_SGS.x, MESH_SGS.y, MESH_SGS.z);
% %% vortex particle parameters
% PART.hp = PART.h_cutoff * MESH_SGS.dx; % a smoothing radius (i.e., a cutoff length or core size)
% PART.hp = PART.h_cutoff * MESH_SGS.dx; % a smoothing radius (i.e., a cutoff length or core size)
% 
% switch SIM.example
%     case 'VortexRings'
%         % particle positions, vorticity, and strengths.  
%         % NOTE: is this initial field divergenge free (del dot omega)? should check this...
%         [xp, wp, ap, PART.nPart] = init_vortexRing(RING, PART, MESH_SGS, ENV);
%                 
%     case 'Turbine'
%         % this type of simulation starts empty, and relies on "emitter" particles
%         xp         = [];
%         wp         = [];
%         ap         = [];
%         PART.nPart = 0;
% 
%     otherwise
%         
% end