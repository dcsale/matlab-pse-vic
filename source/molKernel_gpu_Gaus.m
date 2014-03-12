function f = molKernel_gpu_Gaus(rx, ry, rz, core)

% user parameter
tol = 1e-6;

r_mag  = sqrt( rx^2 + ry^2 + rz^2 );
r_mag2 = r_mag * r_mag;
core2  = core * core;
core3  = core * core * core;
f = 0;
if r_mag > tol
    f = exp(-r_mag2/(2*core2)) / ((2*pi*core3)^(3/2));
end