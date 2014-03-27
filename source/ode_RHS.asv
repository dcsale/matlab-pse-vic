function dx_dt = ode_RHS(t, x, CTRL, SIM, MESH, PART, ENV)

%% bookkeeping for ODE solver
% reshape the array into more meaningful variable names
nVars = 2;
tmp   = reshape(       x, SIM.dim*PART.nPart, nVars);  
xp    = reshape(tmp(:,1),         PART.nPart, SIM.dim)';    % the particle positions
wp    = reshape(tmp(:,2),         PART.nPart, SIM.dim)';    % the particle weights
% write diagnostics     
fprintf(1, '[ode_RHS.m] time = %g\n', t);

%% Vortex-in-Cell algorithm
volp       = PART.hp^3;                             % particle volumes approx equal to the core size cubed

wf         = interp_P2M(SIM, MESH, xp, wp);       	% init a new vorticity field by interpolation from particles (P2M)
    wf_2 = vort_P2M(SIM, MESH, PART.nPart, xp, wp.*volp, PART.hp, 'GPU-v1');    % for comparison, compute the vorticity field by a full Biot-Savart law
uf         = PoissonSolve3D(SIM, MESH, wf);      	% solve Poisson eqn for velocity
uf         = add_freestream(uf, ENV);               % add the freestream velocity (Helmholtz decomposition)
% (more to come here)                           	% compute strain rate for use in diagnostics and/or LES model
% (more to come here)                               % perform some diagnostics on velocity and vorticity fields
up         = interp_M2P(SIM, MESH, xp, uf);       	% interpolate velocity field to particles (M2P) - the Poisson solver should already have extrapolated velocity into the ghost layer
dwf        = vorticity_rhs(SIM, MESH, wf, uf);   	% compute RHS of vorticity transport eqn - this includes vortex stretching and diffusion calculated on mesh (finite differences) - the mesh and fields should be ghosted/extrpolated before calling
dwp        = interp_M2P(SIM, MESH, xp, dwf);      	% interpolate vorticity RHS from mesh to particles
dap_dt_RHS = volp.*dwp;                             % update the particle weights (inviscid)

dx_dt = [reshape(        up', SIM.dim*PART.nPart, 1); 
         reshape(dap_dt_RHS', SIM.dim*PART.nPart, 1)];

end % function