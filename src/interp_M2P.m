function wp = interp_M2P(PART, MESH, xp, wf)
%% =======================================================================%
% assign particle strengths through the M2P interpolation
% ========================================================================%

% tic

% some temp variables
dim = numel(MESH.dx);

dx = MESH.dx(1);
dy = MESH.dx(2);
dz = MESH.dx(3);

xp_x = xp(1,:);
xp_y = xp(2,:);
xp_z = xp(3,:);

xf = MESH.x{1};
yf = MESH.x{2};
zf = MESH.x{3};


wp_x = zeros(1, PART.nPart);
wp_y = zeros(1, PART.nPart);
wp_z = zeros(1, PART.nPart);

wf_x = wf{1};
wf_y = wf{2};
wf_z = wf{3};

isupport = 2;

parfor p = 1:PART.nPart  
% for p = 1:PART.nPart  
        
            %% loop over all the mesh points only within support
            sx = abs( xp_x(p) - xf ) ./ dx; 
            sy = abs( xp_y(p) - yf ) ./ dy; 
            sz = abs( xp_z(p) - zf ) ./ dz;
%             s  = [sx; sy; sz];
            % the particle is in the support of these mesh points
%             xf_s = xf(sx <= isupport);
%             yf_s = yf(sy <= isupport);
%             zf_s = zf(sz <= isupport);
            % the particle is in the support of these mesh points indicies
            ai = find(sx <= isupport, 1, 'first');
            bi = find(sx <= isupport, 1, 'last');
            aj = find(sy <= isupport, 1, 'first');
            bj = find(sy <= isupport, 1, 'last');
            ak = find(sz <= isupport, 1, 'first');           
            bk = find(sz <= isupport, 1, 'last');
            
            wp_sum = zeros(dim, 1);
            for i = ai:bi
                for j = aj:bj
                    for k = ak:bk 
                        W      = interpKernel([sx(i); ...
                                               sy(j); ...
                                               sz(k)]);
%                         wp_sum = wp_sum + W .* [wf{1}(i,j,k); wf{2}(i,j,k); wf{3}(i,j,k)];
                        wp_sum = wp_sum + W .* [wf_x(i,j,k); wf_y(i,j,k); wf_z(i,j,k)];
                    end
                end
            end

wp_x(p) = wp_sum(1);
wp_y(p) = wp_sum(2);
wp_z(p) = wp_sum(3);

% fprintf(1,'[interp_M2P.m] particle %g of %g. \n', p, PART.nPart);
end

wp = [wp_x; wp_y; wp_z];

% user = memory;
% mem  = user.MemUsedMATLAB;
% MESH.NX
% PART.nPart
% toc 
% mem/1e9
% memory

end % function interp_M2P

function W = interpKernel(s)
% r is DIMENSION x 1 array
% W is SCALAR
dim = numel(s);
W = zeros(dim, 1);
for n = 1:dim
    W(n) = kernel_Mp4(s(n));
end
W = prod(W); % reduce to scalar;
end % interpKernel

function W = kernel_Mp4(s)
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