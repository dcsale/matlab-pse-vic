function dwf = vorticity_rhs(SIM, MESH, wf, uf)

% compute the strain and rotation tensors
[S, ~] = tensors(SIM, MESH, wf, uf);

% compute the vortex stretching term (need to extrapolate into ghost layer?)
w_du_x = wf{1}.*S{1,1} + wf{2}.*S{1,2} + wf{3}.*S{1,3};
w_du_y = wf{1}.*S{2,1} + wf{2}.*S{2,2} + wf{3}.*S{2,3};
w_du_z = wf{1}.*S{3,1} + wf{2}.*S{3,2} + wf{3}.*S{3,3};
dwf    = {w_du_x; w_du_y; w_du_z};

% compute the diffusion term (coming soon)

end % function

function [S, R] = tensors(SIM, MESH, wf, uf)

% gradients (tensors) of the velocity and vorticity fields
g_uf = grad_field(SIM, MESH, uf);
g_wf = grad_field(SIM, MESH, wf);

% rate of strain tensor (symmetric)
S = {               g_uf.dudx, (g_uf.dudy+g_uf.dvdx)./2, (g_uf.dudz+g_uf.dwdx)./2; ...
     (g_uf.dvdx+g_uf.dudy)./2,                g_uf.dvdy, (g_uf.dvdz+g_uf.dwdy)./2; ...
     (g_uf.dwdx+g_uf.dudz)./2, (g_uf.dwdy+g_uf.dvdz)./2,                g_uf.dwdz};
 
% rotation 'vorticity' tensor (anti-symmetric)
R = {        0,  wf{3}./2, -wf{2}./2; ...
     -wf{3}./2,         0,  wf{1}./2; ...
      wf{2}./2, -wf{1}./2,         0};


% divergence of velocity and vorticity fields
div_uf = g_uf.dudx + g_uf.dvdy + g_uf.dwdz;
div_wf = g_wf.dudx + g_wf.dvdy + g_wf.dwdz;
% extrapolate the fields into the ghost layers
div_uf = extrapolate_field_scalar(SIM, MESH, div_uf);
div_wf = extrapolate_field_scalar(SIM, MESH, div_wf);
if SIM.DEBUG_LVL == 999
   % plot the divergence fields for debugging purposes
   plot_scalar(div_uf, MESH, 'velocity divergence');
   plot_scalar(div_wf, MESH, 'vorticity divergence');  
end

end % function

function grad = grad_field(SIM, MESH, uf)

%%------------------------------------------------------------------
% O(4) Finite Difference
%-------------------------------------------------------------------

% initialize, accounting for size of ghost layer
grad.dudx = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dvdx = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dwdx = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dudy = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dvdy = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dwdy = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dudz = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dvdz = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
grad.dwdz = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);

% using ghost layers, only central difference formulas are needed.
facx = 1.0/(MESH.dx(1)*12.0);
facy = 1.0/(MESH.dx(2)*12.0);
facz = 1.0/(MESH.dx(3)*12.0);
for k = 1+SIM.mbc:MESH.NX(3)
    for j = 1+SIM.mbc:MESH.NX(2)
        for i = 1+SIM.mbc:MESH.NX(1)
            grad.dudx(i,j,k) = -    facx*uf{1}(i+2,j  ,k  ) ...
                               +8.0*facx*uf{1}(i+1,j  ,k  ) ...
                               -8.0*facx*uf{1}(i-1,j  ,k  ) ...
                               +    facx*uf{1}(i-2,j  ,k  );
         
            grad.dvdx(i,j,k) = -    facx*uf{2}(i+2,j  ,k  ) ...
                               +8.0*facx*uf{2}(i+1,j  ,k  ) ...
                               -8.0*facx*uf{2}(i-1,j  ,k  ) ...
                               +    facx*uf{2}(i-2,j  ,k  );
      
            grad.dwdx(i,j,k) = -    facx*uf{3}(i+2,j  ,k  ) ...
                               +8.0*facx*uf{3}(i+1,j  ,k  ) ...
                               -8.0*facx*uf{3}(i-1,j  ,k  ) ...
                               +    facx*uf{3}(i-2,j  ,k  );
                
            grad.dudy(i,j,k) = -    facy*uf{1}(i  ,j+2,k  ) ...
                               +8.0*facy*uf{1}(i  ,j+1,k  ) ...
                               -8.0*facy*uf{1}(i  ,j-1,k  ) ...
                               +    facy*uf{1}(i  ,j-2,k  );
                
            grad.dvdy(i,j,k) = -    facy*uf{2}(i  ,j+2,k  ) ...
                               +8.0*facy*uf{2}(i  ,j+1,k  ) ...
                               -8.0*facy*uf{2}(i  ,j-1,k  ) ...
                               +    facy*uf{2}(i  ,j-2,k  );
                
            grad.dwdy(i,j,k) = -    facy*uf{3}(i  ,j+2,k  ) ...
                               +8.0*facy*uf{3}(i  ,j+1,k  ) ...
                               -8.0*facy*uf{3}(i  ,j-1,k  ) ...
                               +    facy*uf{3}(i  ,j-2,k  );
                
            grad.dudz(i,j,k) = -    facz*uf{1}(i  ,j  ,k+2) ...
                               +8.0*facz*uf{1}(i  ,j  ,k+1) ...
                               -8.0*facz*uf{1}(i  ,j  ,k-1) ...
                               +    facz*uf{1}(i  ,j  ,k-2);
                
            grad.dvdz(i,j,k) = -    facz*uf{2}(i  ,j  ,k+2) ...
                               +8.0*facz*uf{2}(i  ,j  ,k+1) ...
                               -8.0*facz*uf{2}(i  ,j  ,k-1) ...
                               +    facz*uf{2}(i  ,j  ,k-2);
                
            grad.dwdz(i,j,k) = -    facz*uf{3}(i  ,j  ,k+2) ...
                               +8.0*facz*uf{3}(i  ,j  ,k+1) ...
                               -8.0*facz*uf{3}(i  ,j  ,k-1) ...
                               +    facz*uf{3}(i  ,j  ,k-2);
        end
    end
end
  
end % function
