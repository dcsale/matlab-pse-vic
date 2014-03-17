function dx_dt = ode_RHS(t, x, CTRL, SIM, MESH, PART, ENV)

% reshape the array into more meaningful variable names
nVars = 2;
tmp   = reshape(       x, SIM.dim*PART.nPart, nVars);  
xp    = reshape(tmp(:,1),         PART.nPart, SIM.dim)';    % the particle positions
wp    = reshape(tmp(:,2),         PART.nPart, SIM.dim)';    % the particle weights




% particle volumes
volp  = PART.hp^3; %  approx equal to the core size cubed

%% calculate RHS of vorticity transport eqn using Vortex-in-Cell algorithm

wf         = interp_P2M(MESH, xp, wp);          % init a new vorticity field by interpolation from particles (P2M)
uf         = PoissonSolve3D(SIM, MESH, wf);     % solve Poisson eqn for velocity
% compute strain rate for use in diagnostics and/or LES model
up         = interp_M2P(MESH, xp, uf);          % interpolate velocity field to particles (M2P)
wf_str     = vortex_stretch(SIM, MESH, wf, uf); % compute vortex stretching on mesh (finite differences)
wp_str     = interp_M2P(MESH, xp, wf_str);      % interpolate vortex stretching from mesh to the particles
% wf_diff    = 0;                                 % compute diffusion on mesh (finite differences) - look into the del2 to compute the discrete Laplacian
% wp_diff    = interp_M2P(MESH, xp, wf_diff);   % interpolate vorticity diffusion from mesh to the particles
% dap_dt_RHS = volp.*(wp + wp_str + wp_diff);     % update the particle weights (viscous)
% dap_dt_RHS = volp.*(wp + wp_str);               % update the particle weights (inviscid)
dap_dt_RHS = volp.*wp_str;               % update the particle weights (inviscid)

dx_dt = [reshape(        up', SIM.dim*PART.nPart, 1); 
         reshape(dap_dt_RHS', SIM.dim*PART.nPart, 1)];

%% write diagnostics     
fprintf(1, '[ode_RHS.m] time = %g\n', t);

end % function