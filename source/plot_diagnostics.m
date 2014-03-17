function plot_diagnostics(SIM, MESH, PART, time, cycle, xp, wp, wf, uf, up)
%Make some plots, showing diagnostics of the vortex method

%% plot the fields and particles for basic visualization
if ~isempty(wf)
    plot_field(wf, MESH, 'Vorticity Field')
end
if ~isempty(uf)
    plot_field(uf, MESH, 'Velocity Field')
end
if ~isempty(xp)
    plot_particles(xp, [], MESH, 'Particle Positions')
end
if ~isempty(wp) && ~isempty(xp)
    plot_particles(xp, wp, MESH, 'Particle Vorticity')
end
if ~isempty(up) && ~isempty(xp)
    plot_particles(xp, up, MESH, 'Particle Velocity')
end

%% =======================================================================%
% plot the divergence of the fields
% ========================================================================%



%% =======================================================================%
% Compare the initial vorticity field to the remeshed particle approximation. 
% The interpolated particle representation of field to original field and compare the moments (1st three moments conserved for M'4 kernel)
% ========================================================================%
% % recompute the field with interpP2M operation
% wf_p2m  = interp_P2M(MESH, xp, wp, 'timing');
% dwf_x   = abs( wf{1} - wf_p2m{1} );
% dwf_y   = abs( wf{2} - wf_p2m{2} );
% dwf_z   = abs( wf{3} - wf_p2m{3} );
% dwf_p2m = {dwf_x; dwf_y; dwf_z};
% plot_field(dwf_p2m, MESH, 'Vorticity Field P2M')

% % recompute the field with full Biot-Savart law
% .
% .
% .
% plot_field(dwf, MESH, 'Vorticity Field P2M')

%% =======================================================================%
% compare curl of velocity field to initial vorticity field
% ========================================================================%
%         % compare the original vorticity field to the curl of velocity
%         wf_curl                              = cell(1, 3);   % initialize vector field in 3D
%         % curl requires MESHGRID format, so convert before calling
%         [wf_curl{1}, wf_curl{2}, wf_curl{3}] = curl(permute(MESH.xf, [2,1,3]), ...
%                                                     permute(MESH.yf, [2,1,3]), ...
%                                                     permute(MESH.zf, [2,1,3]), ...
%                                                     permute(vel{1}, [2,1,3]), ...
%                                                     permute(vel{2}, [2,1,3]), ...
%                                                     permute(vel{3}, [2,1,3]));
%         % convert field from MESHGRID to NDGRID format
%         wf_curl{1} = permute(wf_curl{1}, [2,1,3]);
%         wf_curl{2} = permute(wf_curl{2}, [2,1,3]);
%         wf_curl{3} = permute(wf_curl{3}, [2,1,3]);
%         
%         plot_field(wf_curl, MESH, 'Vorticity Field')
%         plot_field(     wf, MESH, 'Vorticity Field')



end

