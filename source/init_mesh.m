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

% now add the ghost layer
Mesh.x{1} = linspace(Mesh.xmin(1)-Mesh.dx(1)*SIM.mbc, Mesh.xmax(1)+Mesh.dx(1)*SIM.mbc, Mesh.NX(1)+2*SIM.mbc);
Mesh.x{2} = linspace(Mesh.xmin(2)-Mesh.dx(2)*SIM.mbc, Mesh.xmax(2)+Mesh.dx(2)*SIM.mbc, Mesh.NX(2)+2*SIM.mbc);
Mesh.x{3} = linspace(Mesh.xmin(3)-Mesh.dx(3)*SIM.mbc, Mesh.xmax(3)+Mesh.dx(3)*SIM.mbc, Mesh.NX(3)+2*SIM.mbc);

% mesh field
% [Mesh.x, Mesh.y, Mesh.z] = ndgrid(x{1}, x{2}, x{3});
[Mesh.xf, Mesh.yf, Mesh.zf] = ndgrid(Mesh.x{1}, Mesh.x{2}, Mesh.x{3});

% Centred domain
centre        = (Mesh.xmin + Mesh.xmax)/2;
Mesh.xf_cen{1} = Mesh.xf - centre(1);
Mesh.xf_cen{2} = Mesh.yf - centre(2);
Mesh.xf_cen{3} = Mesh.zf - centre(3);

%% should the extended domain be initialized here, or every time within the Poisson solver? ... probably here if freespace BC are chosen
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
