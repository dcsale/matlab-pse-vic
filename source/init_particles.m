function [xp, wp, PART] = init_particles(CTRL, SIM, MESH, ENV)

% to-do
% -rotate the vorticity field according to the ring axis inputs
% -add perturbations to the ring according to user inputs
% -for adaptive meshes, this function also needs to setup and return the
% mesh data--but how to figure out the mesh extents before hand? perhaps
% just grow "r" until threshold values are met--use while loop

gamma  = CTRL.Re * ENV.kin_visc;  	% circulation of vortex ring

% % first need to find the zero level set of the vorticity field
% dr = 
% r  = 0;
% while wf_mag > cutoff
%    r = r + dr;  
% end
% % the vorticity field for multiple vortex rings
% 
% for n = 1:nRings
%     % an array of function handles
%      r           = @(x,y,z) (z - CTRL.center_z(n)).^2 + ( sqrt((x - CTRL.center_x(n)).^2 + (y - CTRL.center_y(n)).^2) - CTRL.Rmajor(n) ).^2;
%      wf_mag      = @(r) gamma(n) / (pi * CTRL.Rminor(n)^2) .* exp(-(r./CTRL.Rminor(n)).^2);
%      wf_mag_norm = wf_mag ./ max(max(max(wf_mag)));
%     
% end

% MESH.xmin = [0, 0, 0]; 
% MESH.xmax = [2, 2, 2];
% % mesh coordinates
% MESH.x{1} = linspace(MESH.xmin(1), MESH.xmax(1), MESH.NX(1));
% MESH.y{2} = linspace(MESH.xmin(2), MESH.xmax(2), MESH.NX(2));
% MESH.z{3} = linspace(MESH.xmin(3), MESH.xmax(3), MESH.NX(3));
% % mesh spacing
% MESH.dx(1) = MESH.x{1}(2) - MESH.x{1}(1);
% MESH.dx(2) = MESH.x{2}(2) - MESH.x{2}(1);
% MESH.dx(3) = MESH.x{3}(2) - MESH.x{3}(1);
% % mesh field
% [Mesh.xf, Mesh.yf, Mesh.zf] = ndgrid(MESH.x{1}, MESH.x{2}, MESH.x{3});
    
xp     = [];
wp     = [];
% ap     = [];
nRings = numel(CTRL.Re);
for n = 1:nRings
  
    %% initialize particles at regular grid locations
%     f = (MESH.Z - CTRL.center_z(n)).^2 + ( sqrt((MESH.X - CTRL.center_x(n)).^2 + (MESH.Y - CTRL.center_y(n)).^2) - CTRL.Rmajor(n) ).^2;
%     f = (MESH.Z).^2 + ( sqrt((MESH.X).^2 + (MESH.Y).^2) - CTRL.Rmajor(n) ).^2;
    f = (MESH.xf{3}).^2 + ( sqrt((MESH.xf{1}).^2 + (MESH.xf{2}).^2) - CTRL.Rmajor(n) ).^2;
    wf_mag_init = gamma(n) / (pi * CTRL.Rminor(n)^2) .* exp(-(f./CTRL.Rminor(n)).^2);
    wf_mag_init_norm = wf_mag_init ./ max(max(max(wf_mag_init)));

%     while wf_mag_init_norm(i,j,k) >= cutoff
%         % how do we find the "0 isosurface" in the quickest manner?
%         % more like, this entire function is bollocks.
%         % but seriously, how to define the initial adaptive mesh?
%     end
    
%     createPart = false(numel(MESH.x), numel(MESH.y), numel(MESH.z));
    createPart = false(numel(MESH.x{1}), numel(MESH.x{2}), numel(MESH.x{3}));
    for i = 1:numel(MESH.x{1})
        for j = 1:numel(MESH.x{2})
            for k = 1:numel(MESH.x{3})
                if wf_mag_init_norm(i,j,k) >= SIM.cutoff        % create a vortex particle here
                    createPart(i,j,k) = true;
                else                                        % do NOT create a vortex particle here (we are below the threshold value for vorticity)
                    createPart(i,j,k) = false;
                end
            end
        end
    end
    nPart_tmp = sum(createPart(:));
    nonZeroEntries = find(createPart);
    [ix, iy, iz] = ind2sub([numel(MESH.x{1}) numel(MESH.x{2}) numel(MESH.x{3})], nonZeroEntries);

    xp_tmp = zeros(3, nPart_tmp);
    wp_tmp = zeros(3, nPart_tmp);
    for p = 1:nPart_tmp
        % particle position
%         xp_tmp(1,p) = MESH.X(ix(p), iy(p), iz(p));
%         xp_tmp(2,p) = MESH.Y(ix(p), iy(p), iz(p));
%         xp_tmp(3,p) = MESH.Z(ix(p), iy(p), iz(p));
        xp_tmp(1,p) = MESH.xf{1}(ix(p), iy(p), iz(p));
        xp_tmp(2,p) = MESH.xf{2}(ix(p), iy(p), iz(p));
        xp_tmp(3,p) = MESH.xf{3}(ix(p), iy(p), iz(p));
        % particle vorticity
        theta       = atan2(xp_tmp(2,p), xp_tmp(1,p));
        wp_tmp(1,p) = -wf_mag_init(ix(p), iy(p), iz(p)) * sin(theta);
        wp_tmp(2,p) = +wf_mag_init(ix(p), iy(p), iz(p)) * cos(theta);
        wp_tmp(3,p) = 0;
    end  
%     ap_tmp = wp_tmp .* PART.hp; 

    %% concatenate the vortex rings particles
    trans      = repmat([CTRL.center_x(n); CTRL.center_y(n); CTRL.center_z(n)], [1 nPart_tmp]);
    xp         = [xp (xp_tmp + trans)];
    wp         = [wp wp_tmp.*CTRL.sign(n)];  % this is a quick hack, remove CTRL.sign(n) in future!
%     ap         = [ap ap_tmp.*CTRL.sign(n)];  % this is a quick hack, remove CTRL.sign(n) in future!
    
    PART.nPart = size(xp, 2);
    PART.hp    = SIM.h_cutoff * max(MESH.dx);     % a smoothing radius (i.e., a cutoff length or core size)
end


end % function

