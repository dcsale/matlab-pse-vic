function [fx, fy, fz] = molKernel_gpu_none(rx, ry, rz)

% user parameter
tol = 1e-6;

r_mag = sqrt( rx^2 + ry^2 + rz^2 );
denom = 4*pi*r_mag^3;

fx = 0;
fy = 0;
fz = 0;
if r_mag > tol
    fx = -rx/denom;
    fy = -ry/denom;
    fz = -rz/denom;
end
     
end 

