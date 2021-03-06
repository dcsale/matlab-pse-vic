function wp = interp_M2P(SIM, MESH, xp, wf, varargin)
%% =======================================================================%
% assign particle strengths through the M2P interpolation
% ========================================================================%

if nargin > 4
    do_timing = true;
    tic
else
    do_timing = false;
end

% some temp variables

dx = MESH.dx(1);
dy = MESH.dx(2);
dz = MESH.dx(3);

xp_x = xp(1,:);
xp_y = xp(2,:);
xp_z = xp(3,:);

xf = MESH.x{1};
yf = MESH.x{2};
zf = MESH.x{3};

nPart = size(xp, 2);
wp_x = zeros(1, nPart);
wp_y = zeros(1, nPart);
wp_z = zeros(1, nPart);

wf_x = wf{1};
wf_y = wf{2};
wf_z = wf{3};

isupport = 2; % NOTE: this should not be hard coded, but depends on the selected interpolation kernel

parfor p = 1:nPart  
% for p = 1:nPart  
        
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
            
            wp_sum = zeros(SIM.dim, 1);
            for i = ai:bi
                for j = aj:bj
                    for k = ak:bk 
                        W      = interp_kernel([sx(i); sy(j); sz(k)]);
                        wp_sum = wp_sum + W .* [wf_x(i,j,k); wf_y(i,j,k); wf_z(i,j,k)];
                    end
                end
            end
wp_x(p) = wp_sum(1);
wp_y(p) = wp_sum(2);
wp_z(p) = wp_sum(3);

% fprintf(1,'[interp_M2P.m] particle %g of %g. \n', p, nPart);
end

%% collect the output
wp = [wp_x; wp_y; wp_z];

%% Benchmarking
if do_timing
%     user = memory;
%     mem  = user.MemUsedMATLAB;
    MESH.NX 
    toc 
%     mem/1e9
%     memory
end

