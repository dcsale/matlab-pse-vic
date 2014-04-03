% function [uf_x, uf_y, uf_z] = vel_P2M_reg(mesh, nPart, xp, ap, hp, SIM.runMode_P2M)
function [uf_x, uf_y, uf_z] = vel_P2M_reg(xp, ap, PART, SIM, MESH, ENV)

if isempty(xp) && isempty(ap)
    uf_x = [];
    uf_y = [];
    uf_z = [];
    return
end

if isstruct(MESH) 
    interpType = 'mesh';
    % MESH is the MESH data structure.  Particle induced velocites will be 
    % calculated on the MESH points.
    % mesh spacing
    dx   = MESH.dx; 
    dy   = MESH.dy;
    dz   = MESH.dz;
    Nx   = MESH.Nx;
    Ny   = MESH.Ny;
    Nz   = MESH.Nz;
    xMin = MESH.x(1);
    yMin = MESH.y(1);
    zMin = MESH.z(1);

elseif isnumeric(MESH)
    % MESH is a [N x 3] array containing points in 3D
    % [px1 py1 pz1; px2 py2 pz2; ...] for which the particle 
    % induced velocity will be calculated.
    interpType = 'points';
    nPoints    = size(MESH, 2);

end

tol = 1e-6;



% particle positions
xp_x = xp(1,:);
xp_y = xp(2,:);
xp_z = xp(3,:);
% particle strengths
ap_x = ap(1,:);
ap_y = ap(2,:);
ap_z = ap(3,:);
% particle core sizes
hp2  = PART.hp^2;
hp3  = PART.hp^3;


switch SIM.runMode_P2M
    case 'CPU-v1'
        % Standard serial code.
        
        uf_x  = zeros(Nx, Ny, Nz);
        uf_y  = zeros(Nx, Ny, Nz);
        uf_z  = zeros(Nx, Ny, Nz);
        for i = 1:Nx
            for j = 1:Ny
                for k = 1:Nz
                    rhs_sum = zeros(3, 1); % re-init to clear old values
                    for p = 1:PART.nPart
                        r           = [MESH.x(i) - xp(1,p); 
                                       MESH.y(j) - xp(2,p); 
                                       MESH.z(k) - xp(3,p)];
                        r_mag       = norm(r);
                        if r_mag < tol
                          continue 
                          % particle is too close, set induced velocity to zero
                        end
                        r_mag2      = r_mag^2;
                        r_mag3      = r_mag^3;
                        fac1        = r_mag3 / (4*pi*hp3);
                        fac2        = (r_mag2/hp2 + 5/2) / ((r_mag2/hp2 + 1)^(5/2));
                        fac3        = fac1 * fac2 / r_mag3;
                        rhs         = -1 .* cross(fac3 .* r, ap(:,p));
                        rhs_sum     = rhs_sum + rhs;
                    end
                    uf_x(i,j,k) = rhs_sum(1);
                    uf_y(i,j,k) = rhs_sum(2);
                    uf_z(i,j,k) = rhs_sum(3);
                end
            end
        end
        
    case 'CPU-v2'
        % By clever indexing, we have reduced 3 for loops to a single
        % parfor loop. Here we parallelize the outer loop over the mesh points.
        % Note: a copy of the particle data (xp and ap) is stored on every processor.
                
        uf_x  = zeros(Nx, Ny, Nz);
        uf_y  = zeros(Nx, Ny, Nz);
        uf_z  = zeros(Nx, Ny, Nz);
        % loop over the indices into uf_x, uf_y, uf_z
        parfor idx = 1:(Nx*Ny*Nz);
            [i, j, k] = ind2sub([Nx Ny Nz], idx);
            xf        = xMin + (i-1)*dx;
            yf        = yMin + (j-1)*dy;
            zf        = zMin + (k-1)*dz;
            
            % for each mesh point, sum the contribution from ALL particles
            rhs_sum = zeros(3, 1); % re-init to clear old values
            for p = 1:PART.nPart
                r           = [xf - xp(1,p); 
                               yf - xp(2,p); 
                               zf - xp(3,p)];
                r_mag       = norm(r);
                if r_mag < tol
                  continue 
                  % particle is too close, set induced velocity to zero
                end
                r_mag2      = r_mag^2;
                r_mag3      = r_mag^3;
                fac1        = r_mag3 / (4*pi*hp3);
                fac2        = (r_mag2/hp2 + 5/2) / ((r_mag2/hp2 + 1)^(5/2));
                fac3        = fac1 * fac2 / r_mag3;
                rhs         = -1 .* cross(fac3 .* r, ap(:,p));
                rhs_sum     = rhs_sum + rhs;
            end
            uf_x(idx) = rhs_sum(1);
            uf_y(idx) = rhs_sum(2);
            uf_z(idx) = rhs_sum(3);
        end
        
  case 'GPU-v1'
    
        % transfer the mesh and particle data to the GPU
        dev_hp   = gpuArray( PART.hp );
        dev_xp_x = gpuArray( xp(1, 1:PART.nPart) );
        dev_xp_y = gpuArray( xp(2, 1:PART.nPart) );
        dev_xp_z = gpuArray( xp(3, 1:PART.nPart) );
        dev_ap_x = gpuArray( ap(1, 1:PART.nPart) );
        dev_ap_y = gpuArray( ap(2, 1:PART.nPart) );
        dev_ap_z = gpuArray( ap(3, 1:PART.nPart) );     
        
        switch interpType
            case 'mesh'
                % allocate the field velocity on the GPU
                dev_uf_x = gpuArray.zeros(Nx, Ny, Nz);
                dev_uf_y = gpuArray.zeros(Nx, Ny, Nz);
                dev_uf_z = gpuArray.zeros(Nx, Ny, Nz);

                % loop over the indices into uf_x, uf_y, uf_z
                parfor idx = 1:(Nx*Ny*Nz);
                    [i, j, k] = ind2sub([Nx Ny Nz], idx);
                    xf        = gpuArray( xMin + (i-1)*dx );
                    yf        = gpuArray( yMin + (j-1)*dy );
                    zf        = gpuArray( zMin + (k-1)*dz );

                    % compute all the node to particle distances
                    dev_rx = xf - dev_xp_x;
                    dev_ry = yf - dev_xp_y;
                    dev_rz = zf - dev_xp_z;
%                     % sort the distances
%                     r_mag  = sort( sqrt(dev_rx.^2 + dev_ry.^2 + dev_rz.^2) , 'ascend');
%                     p_idx  = r_mag < 5*PART.hp;
%                     % pass only the "close enough" particles
%                     [dev_uf_x_pc, dev_uf_y_pc, dev_uf_z_pc] = arrayfun(@induce_vel_gpu_v1, dev_rx(p_idx), ...
%                                                                                            dev_ry(p_idx), ...
%                                                                                            dev_rz(p_idx), ...
%                                                                                            dev_hp, ...
%                                                                                            dev_ap_x(p_idx), ...
%                                                                                            dev_ap_y(p_idx), ...
%                                                                                            dev_ap_z(p_idx));
                    [dev_uf_x_pc, dev_uf_y_pc, dev_uf_z_pc] = arrayfun(@induce_vel_gpu_v1, dev_rx, ...
                                                                                           dev_ry, ...
                                                                                           dev_rz, ...
                                                                                           dev_hp, ...
                                                                                           dev_ap_x, ...
                                                                                           dev_ap_y, ...
                                                                                           dev_ap_z);
                    dev_uf_x(idx) = sum(dev_uf_x_pc);
                    dev_uf_y(idx) = sum(dev_uf_y_pc);
                    dev_uf_z(idx) = sum(dev_uf_z_pc);
                end
                    
            case 'points'
                % allocate the field velocity on the GPU
                dev_uf_x = gpuArray.zeros(nPoints, 1);
                dev_uf_y = gpuArray.zeros(nPoints, 1);
                dev_uf_z = gpuArray.zeros(nPoints, 1);
                
                % loop over the indices into uf_x, uf_y, uf_z
                parfor n = 1:nPoints;
                    % compute all the node to particle distances
                    dev_rx = MESH(1,n) - dev_xp_x;
                    dev_ry = MESH(2,n) - dev_xp_y;
                    dev_rz = MESH(3,n) - dev_xp_z;
%                     % sort the distances
%                     r_mag  = sort( sqrt(dev_rx.^2 + dev_ry.^2 + dev_rz.^2) , 'ascend');
%                     p_idx  = r_mag < 5*PART.hp;
%                     
%                     [dev_uf_x_pc, dev_uf_y_pc, dev_uf_z_pc] = arrayfun(@induce_vel_gpu_v1, dev_rx(p_idx), ...
%                                                                                            dev_ry(p_idx), ...
%                                                                                            dev_rz(p_idx), ...
%                                                                                            dev_hp, ...
%                                                                                            dev_ap_x(p_idx), ...
%                                                                                            dev_ap_y(p_idx), ...
%                                                                                            dev_ap_z(p_idx));
                    [dev_uf_x_pc, dev_uf_y_pc, dev_uf_z_pc] = arrayfun(@induce_vel_gpu_v1, dev_rx, ...
                                                                                           dev_ry, ...
                                                                                           dev_rz, ...
                                                                                           dev_hp, ...
                                                                                           dev_ap_x, ...
                                                                                           dev_ap_y, ...
                                                                                           dev_ap_z);
                    dev_uf_x(n) = sum(dev_uf_x_pc);
                    dev_uf_y(n) = sum(dev_uf_y_pc);
                    dev_uf_z(n) = sum(dev_uf_z_pc);
                end
        end
           
        % add the free stream (irrotational) velocity field
        dev_uf_x = dev_uf_x + ENV.velFree;
        dev_uf_y = dev_uf_y + 0;
        dev_uf_z = dev_uf_z + 0;
        
        % transfer the arrays on the GPU back onto the CPU
        uf_x = gather( dev_uf_x );
        uf_y = gather( dev_uf_y );
        uf_z = gather( dev_uf_z );
        
  case 'GPU-v2'
                
  otherwise
    error('[ERROR] Unrecognized input for variable: SIM.runMode');
end

end % function

