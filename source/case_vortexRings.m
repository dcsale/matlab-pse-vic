function case_vortexRings(CTRL, SIM, MESH, ENV)
% simulation stepping
time   = 0;
cycle  = 0;
dt_out = 1/SIM.fps_output;   % interval to write outputs

% initialize the mesh, field, and particles
MESH = init_mesh(SIM, MESH);                % initialize the mesh
wf   = init_field(CTRL, SIM, MESH, ENV);	% init a field on the mesh
[xp_old, wp_old, PART] = remesh_particles(SIM, MESH, wf);   % initialize particle representation of field (step 1 - initialize the particles)

while time < SIM.endtime-dt_out
    
    % is this the correct way to remesh particles, shouldn't xp_old, wp_old be the input arguments!?
    % there is a difference between remeshing from a field, and remeshing
    % from paticles, we want the later in this situation.
%     [xp_old, wp_old, PART] = remesh_particles(SIM, MESH, wf);   % initialize particle representation of field (step 1 - initialize the particles)

      % integrate in time to determine the updated particle positions and updated strengths (particle volumes unchanged since advection is incompressible)
    x_init = [reshape(xp_old', SIM.dim*PART.nPart, 1); 
              reshape(wp_old', SIM.dim*PART.nPart, 1)];
    switch SIM.TSTEPPING
        case 'fixed'
            tspan = time:SIM.dt:max(time+dt_out, time+SIM.dt);          % fixed time stepping
%             Y = ode1(@ode_RHS, tspan, x_init, CTRL, SIM, MESH, PART, ENV);
            Y = ode2(@ode_RHS, tspan, x_init, CTRL, SIM, MESH, PART, ENV);
%             Y = ode4(@ode_RHS, tspan, x_init, CTRL, SIM, MESH, PART, ENV);
        case 'variable'
            tspan = [time, time+dt_out];	% variable time stepping
%             [~, Y] = ode45(@ode_RHS, tspan, x_init, SIM.optionsODE, CTRL, SIM, MESH, PART, ENV);
            [~, Y] = ode113(@ode_RHS, tspan, x_init, SIM.optionsODE, CTRL, SIM, MESH, PART, ENV);          
    end
       
    nVars  = 2;
    tmp    = reshape(Y(end,:)', SIM.dim*PART.nPart, nVars);         % the solution matrix taking the final time, reassign into meaningful variable names
    xp     = reshape( tmp(:,1),         PART.nPart, SIM.dim)';      % the updated particle positions
    wp     = reshape( tmp(:,2),         PART.nPart, SIM.dim)';      % the updated particle weights (vorticity)
    
    %% output & diagnostics (should move this into the RHS function to avoid recomputing all these things)
    % compute the updated fields from the updated particles
    % update the mesh because particles have now moved
    wf = interp_P2M(SIM, MESH, xp, wp);      	% the updated vorticity field interpolated from updated particles (P2M) - this is already calculated in RHS
    
    uf = PoissonSolve3D(SIM, MESH, wf);         % the updated velocity field, solved by Poisson eqn - this is already calculated in RHS
    up = interp_M2P(SIM, MESH, xp, uf);       	% the updated particle velocities, interpolated from velocity field (M2P) - this is already calculated in RHS
    [uf, up] = add_freestream(uf, up, ENV);     % add the freestream velocity (Helmholtz decomposition)
        
    write_outputs(SIM, MESH, PART, time, cycle, xp, wp, wf, uf, up)
%     plot_diagnostics(SIM, MESH, PART, time, cycle, xp, wp, wf, uf, up);
    plot_particles(xp, up, MESH, 'Particle Velocity')
            
    %% prepare for next cycle
    % update the time span for next integration cycle
    cycle = cycle + 1;
    time  = tspan(end);
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
    save_VectorAndMagnitudeVTK_binary(uf{1}, uf{2}, uf{3}, MESH.xf{1}, MESH.xf{2}, MESH.xf{3}, [SIM.outputDir filesep 'Velocity_Field_'  num2str(cycle,'%04.0f') '.vtk'],'vel_field')
end
if SIM.writeVorticityField
    save_VectorAndMagnitudeVTK_binary(wf{1}, wf{2}, wf{3}, MESH.xf{1}, MESH.xf{2}, MESH.xf{3}, [SIM.outputDir filesep 'Vorticity_Field_' num2str(cycle,'%04.0f') '.vtk'],'vort_field')
end 

fprintf(1, 'WRITE DIAGNOSTICS time     = %g \n', time);
fprintf(1, '                  vol_mesh = %g \n', prod(MESH.dx) * prod(MESH.NX));
fprintf(1, '                  vol_part = %g \n', PART.hp^3 * sum(PART.nPart));
fprintf(1, '                  nPart    = %g \n', PART.nPart);
fprintf(1, '                  uf_max   = %g \n', max(max(max( sqrt(   uf{1}.^2 +   uf{2}.^2 +   uf{3}.^2 ) ))) );
fprintf(1, '                  up_max   = %g \n', max(max(max( sqrt( up(1,:).^2 + up(2,:).^2 + up(3,:).^2 ) ))) );
fprintf(1, '                  wf_max   = %g \n', max(max(max( sqrt(   wf{1}.^2 +   wf{2}.^2 +   wf{3}.^2 ) ))) );
fprintf(1, '                  wp_max   = %g \n', max(max(max( sqrt( wp(1,:).^2 + wp(2,:).^2 + wp(3,:).^2 ) ))) );

end % function
