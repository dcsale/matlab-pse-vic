function W = interp_kernel(s)
% r is DIMENSION x 1 array
% W is SCALAR
dim = numel(s);
W = zeros(dim, 1);
for n = 1:dim
    W(n) = kernel_Mp4(s(n));
end
W = prod(W); % reduce to scalar;
end % interp_kernel

function W = kernel_Mp4(s)
% NOTE: M'4 requires the grid to be Cartesian and uniform
s = abs(s);     % the kernel is symmetric
if s > 2
    W = 0;   
elseif s >= 1 && s <=2
    W = 1/2 * (2 - s)^2 * (1 - s);  
elseif s < 1
    W = 1 - 5/2*s^2 + 3/2*s^3;  
else
    error('[kernel_Mp4.m] ERROR: unknown support')
end
end % kernel_Mp4 