function [w_x_pc, w_y_pc, w_z_pc] = induce_vort_gpu_v1(rx, ...
                                                       ry, ...
                                                       rz, ...
                                                       hp, ...
                                                       ap_x, ...
                                                       ap_y, ...
                                                       ap_z)
tol   = 1e-6;

r_mag  = sqrt( rx^2 + ry^2 + rz^2 );
r_mag2 = r_mag^2;
hp2    = hp^2;
hp3    = hp^3;

f = 0;
if r_mag > tol
    f = exp(-r_mag2/(2*hp2)) / (2*pi*hp3)^(3/2);
end

w_x_pc = ap_x * f;
w_y_pc = ap_y * f;
w_z_pc = ap_z * f;

end % function

