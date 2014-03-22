function case_vortexRings(CTRL, SIM, MESH, ENV)
% simulation stepping
time      = 0;
cycle     = 0;

% set intervals to write outputs
dt_out    = 1/SIM.fps_output;

%% initialize the field and mesh
% [xp, wp, ap, nPart] = init_vortexRing(RING, PART, MESH, ENV)


%% initialize the particles particle from a field
MESH = init_mesh(SIM, MESH);                              % init mesh
% wf   = init_field(CTRL, SIM, MESH);                             % init vorticity field
% [xp_old, wp_old, PART] = remesh_particles(SIM, MESH, wf);       % remesh the field onto new particles (step 1 - initialize the particles)

%% initialize the particles directly and following mesh
[xp_old, wp_old, PART] = init_particles(CTRL, SIM, MESH, ENV);    % init particles directly from restart or function (step 1 - initialize the particles)

while time < SIM.endtime-dt_out
    %% The simulation algorithm thus proceeds as follows:
        % Initialize particles 
        % Do time steps n = 0, . . . , T
            % interpolate ? from the particles onto the mesh
            % solve Poisson eqn for v on the mesh using a fast Poisson solver (e.g. multigrid or FFT)
            % interpolate v from the mesh back to the particles ? v
            % compute vortex stretching ? � ?v on the mesh using, e.g., finite differences
            % interpolate the vortex stretching result from the mesh to the particles as a change of vorticity
            % compute vorticity diffusion ? ?? (isotropic and homogeneous) on
            % the particles using PSE (or RW) (alternatively this can also be computed on the mesh), add this vorticity change to the one from the stretching term, and update the particle circulation. Also add the curl of the body force.
            % move (advect) the particles with velocity v and update their positions, leaving the particle volumes unchanged since advection is incompressible.
            % remesh if needed. Generate particles only at mesh nodes where |field| > 0.
        % end time step  

    % update the mesh from the particle positions
    MESH = updateMesh(SIM, MESH, xp_old, PART);

    % update particle positions: move (advect) the particles with the velocity (particle volumes unchanged since advection is incompressible)
    % integrate in time to determine the updated particle positions and updated strengths
% 	tspan = [time, time+dt_out];                % variable time stepping
    tspan = time:SIM.dt:(time+dt_out);          % fixed time stepping
    x_init = [reshape(xp_old', SIM.dim*PART.nPart, 1); 
              reshape(wp_old', SIM.dim*PART.nPart, 1)];

%     [~, Y] = ode45(@ode_RHS, tspan, x_init, SIM.optionsODE, CTRL, SIM, MESH, PART, ENV);
%     [~, Y] = ode113(@ode_RHS, tspan, x_init, SIM.optionsODE, CTRL, SIM, MESH, PART, ENV);
%     [~, Y] = ode23(@ode_RHS, tspan, x_init, SIM.optionsODE, CTRL, SIM, MESH, PART, ENV);
%     [~, Y] = ode15s(@ode_RHS, tspan, x_init, SIM.optionsODE, CTRL, SIM, MESH, PART, ENV);
%     Y = ode1(@ode_RHS, tspan, x_init, CTRL, SIM, MESH, PART, ENV);
    Y = ode4(@ode_RHS, tspan, x_init, CTRL, SIM, MESH, PART, ENV);
    
    nVars  = 2;
    tmp    = reshape(Y(end,:)', SIM.dim*PART.nPart, nVars);         % the solution matrix taking the final time, reassign into meaningful variable names
    xp     = reshape( tmp(:,1),         PART.nPart, SIM.dim)';      % the updated particle positions
    ap     = reshape( tmp(:,2),         PART.nPart, SIM.dim)';      % the updated particle strenghts (circulation)
    wp     = ap./(PART.hp^3);                                       % the updated particle weights (vorticity)
    
    
    %% output & diagnostics
    % compute the updated fields from the updated particles
    % update the mesh because particles have now moved
    wf = interp_P2M(SIM, MESH, xp, wp);      	% the updated vorticity field interpolated from updated particles (P2M) - this is already calculated in RHS
    uf = PoissonSolve3D(SIM, MESH, wf);  	% the updated velocity field, solved by Poisson eqn - this is already calculated in RHS
    up = interp_M2P(MESH, xp, uf);          % the updated particle velocities, interpolated from velocity field (M2P) - this is already calculated in RHS
    
    write_outputs(SIM, MESH, PART, time, cycle, xp, wp, wf, uf, up)
% 	plot_diagnostics(SIM, MESH, PART, time, cycle, xp, wp, wf, uf, up);
    
    
    
    %% prepare for next cycle
    % update the time span for next integration cycle
    cycle = cycle + 1;
    time  = tspan(2);
    % update the particles
    xp_old = xp;
    wp_old = wp;
    
    
end % time

    

end % function

function write_outputs(SIM, MESH, PART, time, cycle, xp, wp, wf, uf, up)
fprintf(1, 'WRITE OUTPUT FILES time = %g \n', time);
if SIM.writeParticles
    write_Point3D(size(xp, 2), xp, wp, up, cycle, SIM.outputDir, 'all')
end
if SIM.writeVelocityField
    save_VectorAndMagnitudeVTK_binary(uf{1}, uf{2}, uf{3}, MESH.xf, MESH.yf, MESH.zf, [SIM.outputDir filesep 'Velocity_Field_'  num2str(cycle,'%04.0f') '.vtk'],'vel_field')
end
if SIM.writeVorticityField
    save_VectorAndMagnitudeVTK_binary(wf{1}, wf{2}, wf{3}, MESH.xf, MESH.yf, MESH.zf, [SIM.outputDir filesep 'Vorticity_Field_' num2str(cycle,'%04.0f') '.vtk'],'vort_field')
end 

fprintf(1, 'WRITE DIAGNOSTICS time     = %g \n', time);
fprintf(1, '                  vol_mesh = %g \n', prod(MESH.dx) * prod(MESH.NX));
fprintf(1, '                  vol_part = %g \n', PART.hp^3 * sum(PART.nPart));

end % function
