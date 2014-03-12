function wf = interp_M2P(PART, MESH, xp, wp)

% tic

% some temp variables
dim  = numel(MESH.dx);
% x    = MESH.x;
% nx   = MESH.NX;
% dx   = MESH.dx;
% xmin = MESH.xmin;
% xmax = MESH.xmax;

%% =======================================================================%
% assign particle strengths through the M2P interpolation
% ========================================================================%
% wp = zeros(3, PART.nPart);
% wp = distributed.zeros(3, PART.nPart);

dx = MESH.dx(1);
dy = MESH.dx(2);
dz = MESH.dx(3);


xp_x = xp(1,:);
xp_y = xp(2,:);
xp_z = xp(3,:);

wp_x = wp(1,:);
wp_y = wp(2,:);
wp_z = wp(3,:);

xf = MESH.x{1};
yf = MESH.x{2};
zf = MESH.x{3};

wf    = cell(dim, 1);
wf{1} = zeros(Mesh.NX);
wf{2} = zeros(Mesh.NX);
wf{3} = zeros(Mesh.NX);

wp_x = zeros(1, PART.nPart);
wp_y = zeros(1, PART.nPart);
wp_z = zeros(1, PART.nPart);

wf_x = wf{1};
wf_y = wf{2};
wf_z = wf{3};

% r_x = zeros(1, numel(nx(1)));
% r_y = zeros(1, numel(nx(2)));
% r_z = zeros(1, numel(nx(3)));

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
%             r   = [r_x; r_y; r_z];
%             r         = [xp(1,p) - x{1}; 
%                          xp(2,p) - x{2}; 
%                          xp(3,p) - x{3}];
%             r         = abs(r);
%             r(1,:)    = r(1,:) ./ dx(1);
%             r(2,:)    = r(2,:) ./ dx(2);
%             r(3,:)    = r(3,:) ./ dx(3);
%             r_support = r <= 2;
%             ai = find(support(1,:), 1, 'first');
%             aj = find(support(2,:), 1, 'first');
%             ak = find(support(3,:), 1, 'first');
%             bi = find(support(1,:), 1, 'last');
%             bj = find(support(2,:), 1, 'last');
%             bk = find(support(3,:), 1, 'last');
%             wp_sum = [0; 0; 0];
%             for i = ai:bi
%                 for j = aj:bj
%                     for k = ak:bk                        
% %                         xf = xmin(1) + (i-1)*dx(1);
% %                         yf = xmin(2) + (j-1)*dx(2);
% %                         zf = xmin(3) + (k-1)*dx(3);
% % 
% %                         r      = [xp(1,p) - xf; 
% %                                   xp(2,p) - yf; 
% %                                   xp(3,p) - zf];
% %                         r      = abs(r);
% %                         r(1,:) = r(1,:) ./ dx(1);
% %                         r(2,:) = r(2,:) ./ dx(2);
% %                         r(3,:) = r(3,:) ./ dx(3);
%                         W      = interpKernel([s]);
% %                         W      = interpKernel(r);
%                         wp_sum = wp_sum + W .* [wf{1}(i,j,k); wf{2}(i,j,k); wf{3}(i,j,k)];   
%                     end
%                 end
%             end

    
%     wp(1,p) = wp_sum(1);
%     wp(2,p) = wp_sum(2);
%     wp(3,p) = wp_sum(3);
    wp_x(p) = wp_sum(1);
    wp_y(p) = wp_sum(2);
    wp_z(p) = wp_sum(3);
    
%     fprintf(1,'[interp_M2P.m] particle %g of %g. \n', p, PART.nPart);

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