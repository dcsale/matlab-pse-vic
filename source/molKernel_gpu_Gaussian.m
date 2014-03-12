function [fx, fy, fz] = molKernel_gpu_Gaussian(rx, ry, rz, hp)

% user parameter
tol = 1e-6;

r_mag  = sqrt( rx^2 + ry^2 + rz^2 );
r_mag3 = r_mag^3;
denom  = 4*pi*r_mag^3;
core3  = hp^3;

fx = 0;
fy = 0;
fz = 0;
if r_mag > tol
    fx = (-rx/denom) * (1 + exp(-r_mag3/core3) - 4*exp(-2*r_mag3/core3));
    fy = (-ry/denom) * (1 + exp(-r_mag3/core3) - 4*exp(-2*r_mag3/core3));
    fz = (-rz/denom) * (1 + exp(-r_mag3/core3) - 4*exp(-2*r_mag3/core3));
end
     
end 
