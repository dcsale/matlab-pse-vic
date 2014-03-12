function [uf_x, uf_y, uf_z] = vel_P2M_sin(mesh, nPart, xp, ap, hp, runMode_P2M)

tol = 1e-6;

% mesh spacing
dx   = mesh.dx; 
dy   = mesh.dy;
dz   = mesh.dz;
Nx   = mesh.Nx;
Ny   = mesh.Ny;
Nz   = mesh.Nz;
xMin = mesh.x(1);
yMin = mesh.y(1);
zMin = mesh.z(1);
% particle positions
xp_x = xp(1,:);
xp_y = xp(2,:);
xp_z = xp(3,:);
% particle strengths
ap_x = ap(1,:);
ap_y = ap(2,:);
ap_z = ap(3,:);
% particle core sizez
hp2  = hp^2;
hp3  = hp^3;


switch runMode_P2M
    case 'CPU-v1'
        % Standard serial code.
        
        uf_x  = zeros(Nx, Ny, Nz);
        uf_y  = zeros(Nx, Ny, Nz);
        uf_z  = zeros(Nx, Ny, Nz);
        for i = 1:Nx
            for j = 1:Ny
                for k = 1:Nz
                    rhs_sum = zeros(3, 1); % re-init to clear old values
                    for p = 1:nPart
                        r           = [mesh.x(i) - xp(1,p); 
                                       mesh.y(j) - xp(2,p); 
                                       mesh.z(k) - xp(3,p)];
                        r_mag       = norm(r);
                        if r_mag < tol
                          continue 
                          % particle is too close, set induced velocity to zero
                        end
                        r_mag3      = r_mag^3;
                        fac1        = 1 ./ (4*pi*r_mag3);
                        rhs         = -1 .* cross(fac1 .* r, ap(:,p));
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
            for p = 1:nPart
                r           = [xf - xp(1,p); 
                               yf - xp(2,p); 
                               zf - xp(3,p)];
                r_mag       = norm(r);
                if r_mag < tol
                  continue 
                  % particle is too close, set induced velocity to zero
                end
                r_mag3      = r_mag^3;
                fac1        = 1 ./ (4*pi*r_mag3);
                rhs         = -1 .* cross(fac1 .* r, ap(:,p));
                rhs_sum     = rhs_sum + rhs;
            end
            uf_x(idx) = rhs_sum(1);
            uf_y(idx) = rhs_sum(2);
            uf_z(idx) = rhs_sum(3);
        end
        
  case 'CPU-v3'      

  case 'GPU-v1'
        
  case 'GPU-v2'
                
  otherwise
    error('[ERROR] Unrecognized input for variable: runMode');
end

end % function

