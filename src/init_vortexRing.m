function [xp, wp, ap, nPart] = init_vortexRing(RING, PART, MESH, ENV)

% to-do
% -rotate the vorticity field according to the ring axis inputs
% -add perturbations to the ring according to user inputs
% -for adaptive meshes, this function also needs to setup and return the
% mesh data--but how to figure out the mesh extents before hand? perhaps
% just grow "r" until threshold values are met--use while loop

gamma  = RING.Re * ENV.kin_visc;  	% circulation of vortex ring
cutoff = 0.01;                      % set vorticity to zero when the field is less than cutoff percent of the field maximum


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
%      r           = @(x,y,z) (z - RING.center_z(n)).^2 + ( sqrt((x - RING.center_x(n)).^2 + (y - RING.center_y(n)).^2) - RING.Rmajor(n) ).^2;
%      wf_mag      = @(r) gamma(n) / (pi * RING.Rminor(n)^2) .* exp(-(r./RING.Rminor(n)).^2);
%      wf_mag_norm = wf_mag ./ max(max(max(wf_mag)));
%     
% end

xp     = [];
wp     = [];
ap     = [];
nRings = numel(RING.Re);
for n = 1:nRings
  
    %% initialize particles at regular grid locations
%     f = (MESH.Z - RING.center_z(n)).^2 + ( sqrt((MESH.X - RING.center_x(n)).^2 + (MESH.Y - RING.center_y(n)).^2) - RING.Rmajor(n) ).^2;
    f = (MESH.Z).^2 + ( sqrt((MESH.X).^2 + (MESH.Y).^2) - RING.Rmajor(n) ).^2;
    wf_mag_init = gamma(n) / (pi * RING.Rminor(n)^2) .* exp(-(f./RING.Rminor(n)).^2);
    wf_mag_init_norm = wf_mag_init ./ max(max(max(wf_mag_init)));

%     while wf_mag_init_norm(i,j,k) >= cutoff
%         % how do we find the "0 isosurface" in the quickest manner?
%         % more like, this entire function is bollocks.
%         % but seriously, how to define the initial adaptive mesh?
%     end
    
    createPart = false(numel(MESH.x), numel(MESH.y), numel(MESH.z));
    for i = 1:numel(MESH.x)
        for j = 1:numel(MESH.y)
            for k = 1:numel(MESH.z)
                if wf_mag_init_norm(i,j,k) >= cutoff   % create a vortex particle here
                    createPart(i,j,k) = true;
                else                                     % do NOT create a vortex particle here (we are below the threshold value for vorticity)
                    createPart(i,j,k) = false;
                end
            end
        end
    end
    nPart_tmp = sum(createPart(:));
    nonZeroEntries = find(createPart);
    [ix, iy, iz] = ind2sub([numel(MESH.x) numel(MESH.y) numel(MESH.z)], nonZeroEntries);

    xp_tmp = zeros(3, nPart_tmp);
    wp_tmp = zeros(3, nPart_tmp);
    for p = 1:nPart_tmp
        % particle position
        xp_tmp(1,p) = MESH.X(ix(p), iy(p), iz(p));
        xp_tmp(2,p) = MESH.Y(ix(p), iy(p), iz(p));
        xp_tmp(3,p) = MESH.Z(ix(p), iy(p), iz(p));
        % particle vorticity
        theta       = atan2(xp_tmp(2,p), xp_tmp(1,p));
        wp_tmp(1,p) = -wf_mag_init(ix(p), iy(p), iz(p)) * sin(theta);
        wp_tmp(2,p) = +wf_mag_init(ix(p), iy(p), iz(p)) * cos(theta);
        wp_tmp(3,p) = 0;
    end  
    ap_tmp = wp_tmp .* PART.hp; 

    %% concatenate the vortex rings particles
    trans      = repmat([RING.center_x(n); RING.center_y(n); RING.center_z(n)], [1 nPart_tmp]);
    xp         = [xp (xp_tmp + trans)];
    wp         = [wp wp_tmp.*RING.sign(n)];  % this is a quick hack, remove RING.sign(n) in future!
    ap         = [ap ap_tmp.*RING.sign(n)];  % this is a quick hack, remove RING.sign(n) in future!
    nPart = size(xp, 2);
end


end % function

