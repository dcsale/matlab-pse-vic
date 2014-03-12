function [uf_x, uf_y, uf_z] = vel_P2M_gpu_v1(xf, yf, zf, xp, wp, h)

% xp_x = xp(1,:);
% xp_y = xp(2,:);
% xp_z = xp(2,:);
% wp_x = wp(1,:);
% wp_y = wp(2,:);
% wp_z = wp(3,:);

% myFUN2 = @(xp_x, xp_y, xp_z, wp_x, wp_y, wp_z) velSum_gpu(xf, ...
%                                                           yf, ...
%                                                           zf, ...
%                                                           xp_x, ...
%                                                           xp_y, ...
%                                                           xp_z, ...
%                                                           wp_x, ...
%                                                           wp_y, ...
%                                                           wp_z, ...
%                                                           h);
                        
% [a2 b2 c2] = myFUN2(gpuArray(1),gpuArray(2),gpuArray(3),gpuArray(4),gpuArray(5),gpuArray(6));
                                                      
% [u_x_cont u_y_cont u_z_cont] = arrayfun(myFUN2, xp_x, xp_y, xp_z, wp_x, wp_y, wp_z);
% 
% uf_x = sum(u_x_cont);
% uf_y = sum(u_y_cont);
% uf_z = sum(u_z_cont);
uf_x = gpuArray(0);
uf_y = gpuArray(0);
uf_z = gpuArray(0);

end

