% [fieldk, fieldx_cen, fieldx_ext, fieldr_ext, NX, vort] = init_field(NXs, SIM.testcase, plot_InitField);
function MESH = init_mesh(SIM, MESH)

%----------------------------------------------------------------------
% Setup domain (_ext is the extended domain)
%----------------------------------------------------------------------

% mesh coordinates (ghost layer added later)
MESH.x{1} = linspace(MESH.xmin(1), MESH.xmax(1), MESH.NX(1));
MESH.x{2} = linspace(MESH.xmin(2), MESH.xmax(2), MESH.NX(2));
MESH.x{3} = linspace(MESH.xmin(3), MESH.xmax(3), MESH.NX(3));

% mesh spacing
MESH.dx(1) = diff( MESH.x{1}(1:2) );
MESH.dx(2) = diff( MESH.x{2}(1:2) );
MESH.dx(3) = diff( MESH.x{3}(1:2) );

% now add the ghost layer
MESH.x{1} = linspace(MESH.xmin(1)-MESH.dx(1)*SIM.mbc, MESH.xmax(1)+MESH.dx(1)*SIM.mbc, MESH.NX(1)+2*SIM.mbc);
MESH.x{2} = linspace(MESH.xmin(2)-MESH.dx(2)*SIM.mbc, MESH.xmax(2)+MESH.dx(2)*SIM.mbc, MESH.NX(2)+2*SIM.mbc);
MESH.x{3} = linspace(MESH.xmin(3)-MESH.dx(3)*SIM.mbc, MESH.xmax(3)+MESH.dx(3)*SIM.mbc, MESH.NX(3)+2*SIM.mbc);

% mesh field
MESH.xf = cell(SIM.dim, 1);
[MESH.xf{1}, MESH.xf{2}, MESH.xf{3}] = ndgrid(MESH.x{1}, MESH.x{2}, MESH.x{3});

% Centred domain
centre        = (MESH.xmin + MESH.xmax)/2;
MESH.xf_cen{1} = MESH.xf{1} - centre(1);
MESH.xf_cen{2} = MESH.xf{2} - centre(2);
MESH.xf_cen{3} = MESH.xf{3} - centre(3);

%% should the extended domain be initialized here, or every time within the Poisson solver? ... probably here if freespace BC are chosen
% % Extended domain (FFT shifted)
% x_ext{1} = MESH.x{1} - CTRL.xmin(1); 
% x_ext{1} = [x_ext{1} -x_ext{1}(end)-MESH.dx(1) -x_ext{1}(end:-1:2)];
% x_ext{2} = MESH.x{2} - CTRL.xmin(2); 
% x_ext{2} = [x_ext{2} -x_ext{2}(end)-MESH.dx(2) -x_ext{2}(end:-1:2)];
% x_ext{3} = MESH.x{3} - CTRL.xmin(3); 
% x_ext{3} = [x_ext{3} -x_ext{3}(end)-MESH.dx(3) -x_ext{3}(end:-1:2)];
% 
% [MESH.xf_ext{1}, MESH.xf_ext{2}, MESH.xf_ext{3}] = ndgrid(x_ext{1},x_ext{2},x_ext{3});
% MESH.rf_ext = sqrt(MESH.xf_ext{1}.^2 + MESH.xf_ext{2}.^2 + MESH.xf_ext{3}.^2);
% 
% if SIM.solve_vel == 2
%     % Wavenumbers for spectral differentiating
%     ks       = 1/MESH.dx(1);
%     k_ext{1} = 2*pi*linspace(-ks/2,ks/2,2*MESH.NX(1)+1);
%     k_ext{1} = fftshift(k_ext{1}(1:end-1));
% 
%     ks       = 1/MESH.dx(2);
%     k_ext{2} = 2*pi*linspace(-ks/2,ks/2,2*MESH.NX(2)+1);
%     k_ext{2} = fftshift(k_ext{2}(1:end-1));
% 
%     ks       = 1/MESH.dx(3);
%     k_ext{3} = 2*pi*linspace(-ks/2,ks/2,2*MESH.NX(3)+1);
%     k_ext{3} = fftshift(k_ext{3}(1:end-1));
% 
%     [MESH.kf{1}, MESH.kf{2}, MESH.kf{3}] = ndgrid(k_ext{1}, k_ext{2}, k_ext{3});
% end

%     clear MESH.x MESH.y MESH.z
        
end % function
