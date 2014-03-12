function dap_dt_RHS = get_RHS(nPart, xp, ap, hp, kin_visc, runMode_RHS)

% particle positions
xp_x = xp(1,:);
xp_y = xp(2,:);
xp_z = xp(3,:);
% particle strengths
ap_x = ap(1,:);
ap_y = ap(2,:);
ap_z = ap(3,:);
% particle volumes
volp = hp^3; %  approx equal to the core size cubed
hp2  = hp^2;
hp4  = hp^4;

switch runMode_RHS
    case 'CPU-v1'
      % Standard serial code.

      dapx_dt_RHS = zeros(1, nPart);
      dapy_dt_RHS = zeros(1, nPart);
      dapz_dt_RHS = zeros(1, nPart);
      % loop over the particles
      for p = 1:nPart;
          % for each particle, sum the contribution from ALL other particles
          rhs_sum = zeros(3, 1); % re-init to clear old values
          for q = 1:nPart
              r       = [xp_x(p) - xp_x(q); 
                         xp_y(p) - xp_y(q); 
                         xp_z(p) - xp_z(q)];
              r_mag   = norm(r);
              r_mag2  = r_mag^2;
              fac1    = (r_mag2 + 5/2*hp2) ./ ((r_mag + hp2).^(5/2));
              trm1    = cross(fac1.*ap(:,p), ap(:,q));
              fac2    = cross(r, ap(:,q));
              fac3    = dot(ap(:,p), fac2);
              fac4    = 3 * (r_mag2 + 7/2*hp2) / ((r_mag2 + hp2).^(7/2));
              trm2    = (fac4 * fac3) .* r;
              fac5    = 105 * kin_visc * hp4 / ((r_mag2 + hp2)^(9/2));
              fac6    = volp.*ap(:,q) - volp.*ap(:,p);
              trm3    = fac5 .* fac6;
              rhs     = -1/(4*pi) .* (trm1 + trm2 + trm3);
              rhs_sum = rhs_sum + rhs;
          end
          dapx_dt_RHS(p) = rhs_sum(1);
          dapy_dt_RHS(p) = rhs_sum(2);
          dapz_dt_RHS(p) = rhs_sum(3);     
      end  
      
    case 'CPU-v2'
        % Vectorized version of standard serial code (got rid of inner for loop)
        
        dapx_dt_RHS = zeros(1, nPart);
        dapy_dt_RHS = zeros(1, nPart);
        dapz_dt_RHS = zeros(1, nPart);
        % loop over the particles
        for p = 1:nPart;
            % for each particle, sum the contribution from ALL other particles
            rx      = xp_x(p) - xp_x;
            ry      = xp_y(p) - xp_y;
            rz      = xp_z(p) - xp_z;
            r       = [rx; ry; rz];
            r_mag   = sqrt( rx.^2 + ry.^2 + rz.^2 );
            r_mag2  = r_mag.^2;
            fac1    = (r_mag2 + 5/2*hp2) ./ ((r_mag + hp2).^(5/2));
            trm1    = cross(ap(:,p) * fac1, ap);
            fac2    = cross(r, ap);
            fac3    = dot( repmat(ap(:,p), 1, nPart), fac2);
            fac4    = 3 * (r_mag2 + 7/2*hp2) ./ ((r_mag2 + hp2).^(7/2));
            trm2    = repmat((fac4 .* fac3), 3, 1) .* r;
            fac5    = 105 .* kin_visc .* hp4 ./ ((r_mag2 + hp2).^(9/2));
            fac6    = volp.*ap - volp.*repmat(ap(:,p), 1, nPart);
            trm3    = repmat(fac5, 3, 1) .* fac6;
            rhs     = -1/(4*pi) .* (trm1 + trm2 + trm3);
                
            dapx_dt_RHS(p) = sum( rhs(1,:) );
            dapy_dt_RHS(p) = sum( rhs(2,:) );
            dapz_dt_RHS(p) = sum( rhs(3,:) );     
        end
            
  case 'CPU-v3'
        % Added a parfor loop to vectorized version of code
        
        dapx_dt_RHS = zeros(1, nPart);
        dapy_dt_RHS = zeros(1, nPart);
        dapz_dt_RHS = zeros(1, nPart);
        % loop over the particles
        parfor p = 1:nPart;
            % for each particle, sum the contribution from ALL other particles
            rx      = xp_x(p) - xp_x;
            ry      = xp_y(p) - xp_y;
            rz      = xp_z(p) - xp_z;
            r       = [rx; ry; rz];
            r_mag   = sqrt( rx.^2 + ry.^2 + rz.^2 );
            r_mag2  = r_mag.^2;
            fac1    = (r_mag2 + 5/2*hp2) ./ ((r_mag + hp2).^(5/2));
            trm1    = cross(ap(:,p) * fac1, ap);
            fac2    = cross(r, ap);
            fac3    = dot( repmat(ap(:,p), 1, nPart), fac2);
            fac4    = 3 * (r_mag2 + 7/2*hp2) ./ ((r_mag2 + hp2).^(7/2));
            trm2    = repmat((fac4 .* fac3), 3, 1) .* r;
            fac5    = 105 .* kin_visc .* hp4 ./ ((r_mag2 + hp2).^(9/2));
            fac6    = volp.*ap - volp.*repmat(ap(:,p), 1, nPart);
            trm3    = repmat(fac5, 3, 1) .* fac6;
            rhs     = -1/(4*pi) .* (trm1 + trm2 + trm3);
                
            dapx_dt_RHS(p) = sum( rhs(1,:) );
            dapy_dt_RHS(p) = sum( rhs(2,:) );
            dapz_dt_RHS(p) = sum( rhs(3,:) );     
        end
        
  case 'GPU-v1'
        % transfer the mesh and particle data to the GPU
        dev_hp   = gpuArray( hp );
        dev_xp_x = gpuArray( xp(1, 1:nPart) );
        dev_xp_y = gpuArray( xp(2, 1:nPart) );
        dev_xp_z = gpuArray( xp(3, 1:nPart) );
        dev_ap_x = gpuArray( ap(1, 1:nPart) );
        dev_ap_y = gpuArray( ap(2, 1:nPart) );
        dev_ap_z = gpuArray( ap(3, 1:nPart) );     
        
        % allocate the particle velocity on the GPU
        dev_dap_dt_x = gpuArray.zeros(1, nPart);
        dev_dap_dt_y = gpuArray.zeros(1, nPart);
        dev_dap_dt_z = gpuArray.zeros(1, nPart);
        
        % loop over the particles
        for p = 1:nPart;
            
            dev_rx = dev_xp_x(p) - dev_xp_x;
            dev_ry = dev_xp_y(p) - dev_xp_y;
            dev_rz = dev_xp_z(p) - dev_xp_z;
            [dev_dap_dt_RHS_x_pc, dev_dap_dt_RHS_y_pc, dev_dap_dt_RHS_z_pc] = arrayfun(@get_RHS_gpu_v1, dev_rx, ...
                                                                                                        dev_ry, ...
                                                                                                        dev_rz, ...
                                                                                                        dev_hp, ...
                                                                                                        dev_ap_x, ...
                                                                                                        dev_ap_y, ...
                                                                                                        dev_ap_z, ...
                                                                                                        dev_ap_x(p), ...
                                                                                                        dev_ap_y(p), ...
                                                                                                        dev_ap_z(p), ...
                                                                                                        kin_visc);
                                                              
            dev_dap_dt_x(p) = sum(dev_dap_dt_RHS_x_pc);
            dev_dap_dt_y(p) = sum(dev_dap_dt_RHS_y_pc);
            dev_dap_dt_z(p) = sum(dev_dap_dt_RHS_z_pc);
         
        end

        % transfer the arrays on the GPU back onto the CPU
        dapx_dt_RHS = gather( dev_dap_dt_x );
        dapy_dt_RHS = gather( dev_dap_dt_y );
        dapz_dt_RHS = gather( dev_dap_dt_z );
        
  case 'GPU-v2'
    
  otherwise
      error('[ERROR] Unrecognized input for variable: runMode');
        
end % switch

dap_dt_RHS = [dapx_dt_RHS; 
              dapy_dt_RHS;
              dapz_dt_RHS];

end % function