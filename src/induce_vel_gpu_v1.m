function [u_x_pc, u_y_pc, u_z_pc] = induce_vel_gpu_v1(rx, ...
                                                      ry, ...
                                                      rz, ...
                                                      hp, ...
                                                      ap_x, ...
                                                      ap_y, ...
                                                      ap_z)
hp2 = hp^2;

r_mag        = sqrt( rx^2 + ry^2 + rz^2 );
r_mag2       = r_mag*r_mag;
c1           = (r_mag2 + 5/2*hp2) / ((r_mag2 + hp2)^(5/2));
[ux, uy, uz] = crossf(c1*rx, c1*ry, c1*rz, ap_x, ap_y, ap_z);

u_x_pc       = -1/(4*pi) * ux;
u_y_pc       = -1/(4*pi) * uy;
u_z_pc       = -1/(4*pi) * uz;

end % function
