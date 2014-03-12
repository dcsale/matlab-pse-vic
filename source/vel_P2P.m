function up = vel_P2P(xp, ap, PART, SIM, ENV)

% particle positions
xp_x = xp(1,:);
xp_y = xp(2,:);
xp_z = xp(3,:);
% particle strengths
ap_x = ap(1,:);
ap_y = ap(2,:);
ap_z = ap(3,:);
% particle volumes
PART.hp2  = PART.hp.^2;

switch SIM.runMode_P2P
    case 'CPU-v1'
      % Standard serial code.

      up_x = zeros(1, PART.nPart);
      up_y = zeros(1, PART.nPart);
      up_z = zeros(1, PART.nPart);
      % loop over the particles
      for p = 1:PART.nPart;
          % for each particle, sum the contribution from ALL other particles
          vel_sum = zeros(3, 1); % re-init to clear old values
          for q = 1:PART.nPart
              r       = [xp_x(p) - xp_x(q); 
                         xp_y(p) - xp_y(q); 
                         xp_z(p) - xp_z(q)];
              r_mag   = norm(r);
              r_mag2  = r_mag^2;
              c1      = (r_mag2 + 5/2*PART.hp2) ./ ((r_mag2 + PART.hp2).^(5/2));
              vel     = -1/(4*pi) .* cross(c1.*r, ap(:,q));
              vel_sum = vel_sum + vel;
          end
          up_x(p) = vel_sum(1);
          up_y(p) = vel_sum(2);
          up_z(p) = vel_sum(3);     
      end  
      
    case 'CPU-v2'
        % Vectorized version of standard serial code (got rid of inner for loop)
        
        up_x = zeros(1, PART.nPart);
        up_y = zeros(1, PART.nPart);
        up_z = zeros(1, PART.nPart);
        for p = 1:PART.nPart
            rx      = xp_x(p) - xp_x;
            ry      = xp_y(p) - xp_y;
            rz      = xp_z(p) - xp_z;
            r_mag   = sqrt( rx.^2 + ry.^2 + rz.^2 );
            r_mag2  = r_mag.^2;
            c1      = (r_mag2 + 5/2*PART.hp2) ./ ((r_mag2 + PART.hp2).^(5/2));
            vel     = -1/(4*pi) .* cross([c1.*rx; c1.*ry; c1.*rz], ap);
            
            up_x(p) = sum( vel(1,:) );
            up_y(p) = sum( vel(2,:) );
            up_z(p) = sum( vel(3,:) );
        end
      
  case 'CPU-v3'
        % Added parfor loop to vectorized version of code
        
        up_x = zeros(1, PART.nPart);
        up_y = zeros(1, PART.nPart);
        up_z = zeros(1, PART.nPart);
        parfor p = 1:PART.nPart
            rx      = xp_x(p) - xp_x;
            ry      = xp_y(p) - xp_y;
            rz      = xp_z(p) - xp_z;
            r_mag   = sqrt( rx.^2 + ry.^2 + rz.^2 );
            r_mag2  = r_mag.^2;
            c1      = (r_mag2 + 5/2*PART.hp2) ./ ((r_mag2 + PART.hp2).^(5/2));
            vel     = -1/(4*pi) .* cross([c1.*rx; c1.*ry; c1.*rz], ap);
            
            up_x(p) = sum( vel(1,:) );
            up_y(p) = sum( vel(2,:) );
            up_z(p) = sum( vel(3,:) );
        end
        
  case 'GPU-v1'
        % transfer the mesh and particle data to the GPU device
        dev_hp   = gpuArray( PART.hp );
        dev_xp_x = gpuArray( xp(1, 1:PART.nPart) );
        dev_xp_y = gpuArray( xp(2, 1:PART.nPart) );
        dev_xp_z = gpuArray( xp(3, 1:PART.nPart) );
        dev_ap_x = gpuArray( ap(1, 1:PART.nPart) );
        dev_ap_y = gpuArray( ap(2, 1:PART.nPart) );
        dev_ap_z = gpuArray( ap(3, 1:PART.nPart) );     
        
        % allocate the particle velocity on the GPU device
        dev_up_x  = gpuArray.zeros(1, PART.nPart);
        dev_up_y  = gpuArray.zeros(1, PART.nPart);
        dev_up_z  = gpuArray.zeros(1, PART.nPart);
        
        % loop over the particles
        for p = 1:PART.nPart;
            
            dev_rx = dev_xp_x(p) - dev_xp_x;
            dev_ry = dev_xp_y(p) - dev_xp_y;
            dev_rz = dev_xp_z(p) - dev_xp_z;
            [dev_up_x_pc, dev_up_y_pc, dev_up_z_pc] = arrayfun(@induce_vel_gpu_v1, dev_rx, ...
                                                                                   dev_ry, ...
                                                                                   dev_rz, ...
                                                                                   dev_hp, ...
                                                                                   dev_ap_x, ...
                                                                                   dev_ap_y, ...
                                                                                   dev_ap_z);
                                                              
            dev_up_x(p) = sum(dev_up_x_pc);
            dev_up_y(p) = sum(dev_up_y_pc);
            dev_up_z(p) = sum(dev_up_z_pc);
         
        end
        
        % add the free stream (irrotational) velocity field
        dev_up_x = dev_up_x + ENV.velFree(1);
        dev_up_y = dev_up_y + ENV.velFree(2);
        dev_up_z = dev_up_z + ENV.velFree(3);
        
        % transfer the arrays on the GPU back onto the CPU host
        up_x = gather( dev_up_x );
        up_y = gather( dev_up_y );
        up_z = gather( dev_up_z );
        
  case 'GPU-v2'
        % transfer the mesh and particle data to the GPU device
        dev_hp   = gpuArray( PART.hp );
        dev_xp_x = gpuArray( xp(1, 1:PART.nPart) );
        dev_xp_y = gpuArray( xp(2, 1:PART.nPart) );
        dev_xp_z = gpuArray( xp(3, 1:PART.nPart) );
        dev_ap_x = gpuArray( ap(1, 1:PART.nPart) );
        dev_ap_y = gpuArray( ap(2, 1:PART.nPart) );
        dev_ap_z = gpuArray( ap(3, 1:PART.nPart) );     
        
        % allocate the particle velocity on the GPU device
        dev_up_x  = gpuArray.zeros(1, PART.nPart);
        dev_up_y  = gpuArray.zeros(1, PART.nPart);
        dev_up_z  = gpuArray.zeros(1, PART.nPart);
        
        % loop over the particles
        for p = 1:PART.nPart;
            
            dev_rx = dev_xp_x(p) - dev_xp_x;
            dev_ry = dev_xp_y(p) - dev_xp_y;
            dev_rz = dev_xp_z(p) - dev_xp_z;
            % particle contributions to particle velocity
            [dev_up_x_pc, dev_up_y_pc, dev_up_z_pc] = arrayfun(@induce_vel_gpu_v2, dev_rx, ...
                                                                                   dev_ry, ...
                                                                                   dev_rz, ...
                                                                                   dev_hp, ...
                                                                                   dev_ap_x, ...
                                                                                   dev_ap_y, ...
                                                                                   dev_ap_z);
                                                              
            dev_up_x(p) = sum(dev_up_x_pc);
            dev_up_y(p) = sum(dev_up_y_pc);
            dev_up_z(p) = sum(dev_up_z_pc);
         
        end
        
        % add the free stream (irrotational) velocity field
        dev_up_x = dev_up_x + ENV.velFree(1);
        dev_up_y = dev_up_y + ENV.velFree(2);
        dev_up_z = dev_up_z + ENV.velFree(3);
        
        % transfer the arrays on the GPU back onto the CPU host
        up_x = gather( dev_up_x );
        up_y = gather( dev_up_y );
        up_z = gather( dev_up_z );
    
  otherwise
      error('[ERROR] Unrecognized input for variable: runMode');
        
end % switch

up = [up_x; 
      up_y;
      up_z];

end % function
