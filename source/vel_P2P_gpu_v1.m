function [uf_x_pc, uf_y_pc, uf_z_pc] = velSum_gpu(xf, ...
                                                  yf, ...
                                                  zf, ...
                                                  h, ...
                                                  xp_x, ...
                                                  xp_y, ...
                                                  xp_z, ...
                                                  wp_x, ...
                                                  wp_y, ...
                                                  wp_z)
% == Inputs ==
% xf   is scalar
% yf   is scalar
% zf   is scalar
% h    is scalar
% xp_x is size(nPart, 1)
% xp_y is size(nPart, 1)
% xp_z is size(nPart, 1)
% wp_x is size(nPart, 1)
% wp_y is size(nPart, 1)
% wp_z is size(nPart, 1)

% == Outputs ==
% uf_x_pc is size(1, nPart)
% uf_y_pc is size(1, nPart)
% uf_z_pc is size(1, nPart)

% h3      = h.^3;
% nPart   = gpuArray.numel(xp_x);
% r       = gpuArray.zeros(3, 1);
% wp      = gpuArray.zeros(3, 1);
% uf_x_pc = gpuArray.zeros(nPart, 1);
% uf_y_pc = gpuArray.zeros(nPart, 1);
% uf_z_pc = gpuArray.zeros(nPart, 1);
% 
% for n = 1:nPart
%   r(1)       = xf - xp_x(n); 
%   r(2)       = yf - xp_y(n); 
%   r(3)       = zf - xp_z(n); 
%   wp(1)      = wp_x(n); 
%   wp(2)      = wp_y(n); 
%   wp(3)      = wp_z(n);
%   K          = molKernel_gpu_none(r);
%   vel_pc     = crossf(K, wp*h3);
%   uf_x_pc(n) = vel_pc(1);
%   uf_y_pc(n) = vel_pc(2);
%   uf_z_pc(n) = vel_pc(3);
% end

h3      = h.^3;
r       = gpuArray.zeros(3, 1);
% wp      = gpuArray.zeros(3, 1);
% uf_x_pc = gpuArray.zeros(nPart, 1);
% uf_y_pc = gpuArray.zeros(nPart, 1);
% uf_z_pc = gpuArray.zeros(nPart, 1);

  r  = [xf - xp_x; 
        yf - xp_y; 
        zf - xp_z];
  wp = [wp_x; 
        wp_y; 
        wp_z];
  K          = molKernel_gpu_none(r);
  vel_pc     = crossf(K, wp.*h3);
  uf_x_pc = vel_pc(1);
  uf_y_pc = vel_pc(2);
  uf_z_pc = vel_pc(3);


end % function

