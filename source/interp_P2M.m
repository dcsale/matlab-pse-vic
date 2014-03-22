function wf = interp_P2M(SIM, MESH, xp, wp, varargin)
%% =======================================================================%
% assign particle strengths through the P2M interpolation
% ========================================================================%

if nargin > 3
    do_timing = true;
    tic
else
    do_timing = false;
end

% some temp variables
% dim  = numel(MESH.dx);
% xMin = MESH.xmin(1);
% yMin = MESH.xmin(2);
% zMin = MESH.xmin(3);

nPart = size(xp, 2);
xp_x  = xp(1,:);
xp_y  = xp(2,:);
xp_z  = xp(3,:);
wp_x  = wp(1,:);
wp_y  = wp(2,:);
wp_z  = wp(3,:);

xf = MESH.x{1};
yf = MESH.x{2};
zf = MESH.x{3};

% Nx = MESH.NX(1);
% Ny = MESH.NX(2);
% Nz = MESH.NX(3);

dx = MESH.dx(1);
dy = MESH.dx(2);
dz = MESH.dx(3);

% wf    = cell(dim, 1);
% wf{1} = zeros(MESH.NX);
% wf{2} = zeros(MESH.NX);
% wf{3} = zeros(MESH.NX);
% wf_x  = wf{1};
% wf_y  = wf{2};
% wf_z  = wf{3};
% wf_x = zeros(MESH.NX);
% wf_y = zeros(MESH.NX);
% wf_z = zeros(MESH.NX);
wf_x = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
wf_y = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);
wf_z = zeros(MESH.NX(1)+2*SIM.mbc, MESH.NX(2)+2*SIM.mbc, MESH.NX(3)+2*SIM.mbc);

isupport = 2; % NOTE: this should not be hard coded, but depends on the selected interpolation kernel

%% version 1
% 
% parfor p = 1:nPart  
for p = 1:nPart  
        
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
            
            for i = ai:bi
                for j = aj:bj
                    for k = ak:bk 
                        W           = interp_kernel([sx(i); sy(j); sz(k)]);
                        wf_x(i,j,k) = wf_x(i,j,k) + W*wp_x(p);
                        wf_y(i,j,k) = wf_y(i,j,k) + W*wp_y(p);
                        wf_z(i,j,k) = wf_z(i,j,k) + W*wp_z(p);
                    end
                end
            end
    
%     fprintf(1,'[interp_P2M.m] particle %g of %g. \n', p, nPart);

end
 
%% collect the output
wf = {wf_x; wf_y; wf_z};

%% for benchmarking
if do_timing
%     user = memory;
%     mem  = user.MemUsedMATLAB;
    MESH.NX 
    toc 
%     mem/1e9
%     memory
end

end
