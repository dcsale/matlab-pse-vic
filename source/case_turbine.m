function case_turbine(CTRL, SIM, MESH, ENV)

%% user parameters
dt_emit = 0.15;     % emit particles every [seconds] -- it's best to choose this value so that the blade tip does not cross more than 1 cell per time step -- or the FAST developers suggest that you choose about 200 time steps per rotor revolution.
theta   = 0;        % initial value. shaft angular position.


%%

% simulation stepping
time      = 0;
cycle     = 0;

% simulation sub-stepping
emit_t     = time:dt_emit:SIM.endtime;
emit_cycle = 0;

% set intervals to write outputs
dt_out    = 1/SIM.fps_output;

% pre-allocate particles & fields
% xp     = [];
% ap     = [];
xp_old = [];
ap_old = [];

% xp_wake     = [];
% ap_wake     = [];
xp_old_wake = [];
ap_old_wake = [];
% wp_old_wake = [];
% up_old_wake = [];

if CTRL.INFLOW_TURB
    % xp_inflow     = [];
    % ap_inflow     = [];
    xp_old_inflow = [];
    ap_old_inflow = [];
    % wp_old_inflow = [];
    % up_old_inflow = [];
end

Coord         = cell(CTRL.NUM_BLADES, 1);
xp_emit_rotor = cell(CTRL.NUM_BLADES, 1);

while time < SIM.endtime-dt_out

%     tspan = [time, time+dt_out];        % variable time stepping
%     tspan = 0:0.1:SIM.endtime;          % fixed time stepping

    %% compute the rotor speed and blade coordinates
    theta_init  = theta;
    switch SIM.TSTEPPING
        case 'fixed'
            tspan = time:SIM.dt:max(time+dt_out, time+SIM.dt);
            Y     =   ode2(@shaftSpeed, tspan, theta_init, SIM.optionsODE, CTRL.ROT_SPD);
        case 'variable'
            tspan = [time, time+dt_out];
            [~, Y] = ode45(@shaftSpeed, tspan, theta_init, SIM.optionsODE, CTRL.ROT_SPD);
    end
    theta_old   = Y(1); % solution at tspan(1)
    theta       = Y(2); % solution at tspan(end)

    azim_shaft_old  = mod(radtodeg(theta_old), 360);
    azim_shaft      = mod(radtodeg(theta),     360);
    azim_blades_old = mod(azim_shaft_old + 360/CTRL.NUM_BLADES .* (0:CTRL.NUM_BLADES-1)', 360);
    azim_blades     = mod(azim_shaft     + 360/CTRL.NUM_BLADES .* (0:CTRL.NUM_BLADES-1)', 360);
    d_dt_theta_old  = shaftSpeed(tspan(1), theta_old, CTRL.ROT_SPD);
    d_dt_theta      = shaftSpeed(tspan(2),     theta, CTRL.ROT_SPD);

    %% compute the emitted particle positions and strengths
    xi = cosspace(0, 1, CTRL.NUM_SEC, 'both'); % blade element spacing
    for n = 1:CTRL.NUM_BLADES
        % a structure array, containing transformations to each blade
        Coord(n)         = defineCoordSystems(CTRL, azim_blades(n));
        xp_emit_rotor(n) = interparc(xi, ...
                                     [Coord(n).Origin.Blade(1) Coord(n).Tip.Blade(1)], ...
                                     [Coord(n).Origin.Blade(2) Coord(n).Tip.Blade(2)], ...
                                     [Coord(n).Origin.Blade(3) Coord(n).Tip.Blade(3)], ...
                                     'linear');
    end
    % the rotor particle should be updated every time step,
    % regardless of whether particles were emitted or not
    xp_rotor = [xp_emit_rotor{1}', ...
                xp_emit_rotor{2}', ...
                xp_emit_rotor{3}'];
    ap_rotor = [repmat(Coord(1).Tran.BG * [0; 0; 1], 1, CTRL.NUM_SEC), ...
                repmat(Coord(2).Tran.BG * [0; 0; 1], 1, CTRL.NUM_SEC), ...
                repmat(Coord(3).Tran.BG * [0; 0; 1], 1, CTRL.NUM_SEC)];

    %% =======================================================================%
    % create particles subject to vorticity transport equation (VTE)
    % ========================================================================%            
    if tspan(2) > emit_t(emit_cycle+1)
        emit_cycle = emit_cycle + 1;

        % wake
        xp_emit_wake = xp_rotor;
        circulation  = 1;           % give a unit circulation in the blade coordinate system
        ap_emit_wake = [repmat(Coord(1).Tran.BG * [0; 0; circulation], 1, CTRL.NUM_SEC), ...
                        repmat(Coord(2).Tran.BG * [0; 0; circulation], 1, CTRL.NUM_SEC), ...
                        repmat(Coord(3).Tran.BG * [0; 0; circulation], 1, CTRL.NUM_SEC)];

                    
        if CTRL.INFLOW_TURB
            
            
        else
            
            
        end
        
        
        % inflow turbulence
        [my, mz] = meshgrid(linspace(             -CTRL.ROTOR_DIA,               CTRL.ROTOR_DIA, 5), ...
                            linspace(CTRL.HUB_HT - CTRL.ROTOR_DIA, CTRL.HUB_HT + CTRL.ROTOR_DIA, 5));
        xp_emit_inflow = [-CTRL.ROTOR_DIA/2.*ones(1,numel(my)); my(:)'; mz(:)'];
        ap_emit_inflow = rand(3, size(xp_emit_inflow, 2));

        % consolidate particles            
        xp        = [xp_old       , xp_emit_wake  , xp_emit_inflow];
        ap        = [ap_old       , ap_emit_wake  , ap_emit_inflow];
        xp_wake   = [xp_old_wake  , xp_emit_wake             ];
%         ap_wake   = [ap_old_wake  , ap_emit_wake             ];
        xp_inflow = [xp_old_inflow, xp_emit_inflow           ];
%         ap_inflow = [ap_old_inflow, ap_emit_inflow           ];
    else
        % consolidate particles
        xp        = xp_old;
        ap        = ap_old;
        xp_wake   = xp_old_wake;
%         ap_wake   = ap_old_wake;
        xp_inflow = xp_old_inflow;
%         ap_inflow = ap_old_inflow;
    end
    PART.nPart        = size(xp, 2);
    PART.nPart_wake   = size(xp_wake  , 2);
    PART.nPart_inflow = size(xp_inflow, 2);
%             PART.nPart_wake   = size([xp_old_wake  , xp_emit_wake  ], 2);
%             PART.nPart_inflow = size([xp_old_inflow, xp_emit_inflow], 2);




    % pass in ALL the particles (inflow, wake, rotor, image), and the points to interpolate to
    % should just preallocate up_rotor here !!! check preallocation in this file to speedup code
    [up_x_rotor, ...
     up_y_rotor, ...
     up_z_rotor] = vel_P2M_reg(xp, ap, PART, SIM, xp_rotor, ENV);
    up_rotor     = [up_x_rotor'; ...
                    up_y_rotor'; ...
                    up_z_rotor'];

    [up_x_wake, ...
     up_y_wake, ...
     up_z_wake] = vel_P2M_reg(xp, ap, PART, SIM, xp_wake, ENV);
    up_wake     = [up_x_wake'; ...
                   up_y_wake'; ...
                   up_z_wake'];

    [up_x_inflow, ...
     up_y_inflow, ...
     up_z_inflow] = vel_P2M_reg(xp, ap, PART, SIM, xp_inflow, ENV);
    up_inflow     = [up_x_inflow'; ...
                     up_y_inflow'; ...
                     up_z_inflow'];

    [up_x, ...
     up_y, ...
     up_z] = vel_P2M_reg(xp, ap, PART, SIM, xp, ENV);
    up     = [up_x'; ...
              up_y'; ...
              up_z'];

    % now with the velocity at the bound particles, lookup the
    % airfoil data to calculate the circulation strength
    % ???

    % integrate in time to determine the new particle positions and strengths
    x_init  = [reshape(xp',3*PART.nPart,1); 
               reshape(ap',3*PART.nPart,1)];
    [~, Y]  = ode113(@get_RHS_ode, tspan, x_init, SIM.optionsODE, PART, ENV, SIM);
    tmp     = reshape(Y(end,:)', 3*PART.nPart, 2);  % the solution matrix
    xp      = reshape( tmp(:,1),   PART.nPart, 3)';
    ap      = reshape( tmp(:,2),   PART.nPart, 3)';

    % instead use a push, pop, stack idea here -- this indexing becomes clumbsy
    xp_wake   = xp(:, PART.nPart + 1 - PART.nPart_inflow - PART.nPart_wake : PART.nPart  - PART.nPart_inflow);
    ap_wake   = ap(:, PART.nPart + 1 - PART.nPart_inflow - PART.nPart_wake : PART.nPart  - PART.nPart_inflow);
    xp_inflow = xp(:, PART.nPart + 1 - PART.nPart_inflow                        : end);
    ap_inflow = ap(:, PART.nPart + 1 - PART.nPart_inflow                        : end);

    %% =======================================================================%
    % update particle velocities: velocity at particle locations
    % ========================================================================%
%             up = vel_P2P(xp, ap, PART, SIM, ENV);    
    %NOTE: "up" is already computed in the get_RHS_ode function! how can we avoid recomputing it here? we need additional outputs from ode113 above...                         

    % =======================================================================%
    % update the mesh to follow the support of vorticity
    % ========================================================================%
    if MESH.adaptive
        MESH = updateMesh(MESH, CTRL.ROTOR_DIA/2, xp, PART);
    end

    % =======================================================================%
    % update velocity field: velocity at mesh point locations
    % ========================================================================%
%             [uf_x, uf_y, uf_z] = vel_P2M_reg(xp, ap, PART, SIM, MESH, ENV);

    % =======================================================================%
    % update vorticity field: vorticity at mesh point locations
    % ========================================================================%
%             [wf_x, wf_y, wf_z] = vort_P2M(MESH, PART.nPart, xp, ap, PART.hp, SIM.runMode_P2M);

    % =======================================================================%
    % write the output files
    % ========================================================================%
    wp        = ap        ./ PART.hp;
    wp_rotor  = ap_rotor  ./ PART.hp;
    wp_wake   = ap_wake   ./ PART.hp;
    wp_inflow = ap_inflow ./ PART.hp;
    if SIM.writeParticles
        write_Point3D( size(xp_rotor , 2), xp_rotor , wp_rotor , up_rotor , cycle, SIM.outputDir, 'rotor')
        write_Point3D( size(xp_wake  , 2), xp_wake  , wp_wake  , up_wake  , cycle, SIM.outputDir, 'wake')
        write_Point3D( size(xp_inflow, 2), xp_inflow, wp_inflow, up_inflow, cycle, SIM.outputDir, 'inflow')
        write_Point3D( size(xp       , 2), xp       , wp       , up       , cycle, SIM.outputDir, 'all')
%                 write_Point3D( size(xp_wake  , 2), [xp_old_wake  , xp_wake  ], [wp_old_wake  , wp_wake  ], [up_old_wake  , up_wake  ], cycle, SIM.outputDir, 'wake')
%                 write_Point3D( size(xp_inflow, 2), [xp_old_inflow, xp_inflow], [wp_old_inflow, wp_inflow], [up_old_inflow, up_inflow], cycle, SIM.outputDir, 'inflow')
    end
    if SIM.writeVelocityField
        save_VectorAndMagnitudeVTK_binary(uf_x, uf_y, uf_z, MESH.X, MESH.Y, MESH.Z, [SIM.outputDir filesep 'Velocity_Field_'  num2str(cycle,'%04.0f') '.vtk'],'vel_field')
    end
    if SIM.writeVorticityField
        save_VectorAndMagnitudeVTK_binary(wf_x, wf_y, wf_z, MESH.X, MESH.Y, MESH.Z, [SIM.outputDir filesep 'Vorticity_Field_' num2str(cycle,'%04.0f') '.vtk'],'vort_field')
    end

    % update the time span for next integration cycle
    cycle = cycle + 1;
    time  = tspan(2);
    % and store any values for next cycle
    xp_old         = xp;
    ap_old         = ap;
    PART.nPart_old = PART.nPart;

    xp_old_wake         = xp_wake;
    ap_old_wake         = ap_wake;
    wp_old_wake         = wp_wake;
    up_old_wake         = up_wake;
    PART.nPart_old_wake = PART.nPart_wake;

    xp_old_inflow         = xp_inflow;
    ap_old_inflow         = ap_inflow;
    wp_old_inflow         = wp_inflow;
    up_old_inflow         = up_inflow;
    PART.nPart_old_inflow = PART.nPart_inflow;

    % =======================================================================%
    % create some quick visuals for debugging purposes
    % ========================================================================%
    if ~exist('fig1') 
    fig1 = figure('name',  'transients', ...
                  'color', 'white', ...
                  'units', 'normalized');
    end
    if ~exist('fig2')
    fig2 = figure('name',  '3D Visualization', ...
                  'color', 'white', ...
                  'units', 'normalized');
    end
    output_exampleTurbine(SIM, MESH, PART, ENV, Coord, CTRL, ...
                          fig1, fig2, ...                          
                          tspan, ...
                          xp, ...
                          xp_emit_rotor, ...
                          azim_blades_old, ...
                          azim_blades, ...
                          d_dt_theta_old, ...
                          d_dt_theta)
%             makeFigures_v2(MESH, CTRL, xp, ap, up, uf_x, uf_y, uf_z, wf_x, wf_y, wf_z, SIM.DEBUG_LVL)   

end
        
end

function d_dt_theta = shaftSpeed(t,theta,rpm)
%Returns d_dt_theta(t) = shaft angular velocity [radians/second]
rps = rpm * pi/30; % shaft speed [radians/second]

t_trans = 2;
t_start = 1;
t_ss    = 10;
t_up    = t_start + t_trans;
t_down  = t_up + t_ss;

if     t < t_start
    % shut-down
    d_dt_theta = 0;
    
elseif t > t_start && t <= t_up
    % ramp up
%     d_dt_theta = rps/t_trans * t;
    d_dt_theta = max(+(rps/t_trans)*(t - t_up) + rps, 0);
    
elseif t > t_up && t <= t_down
    % steady state
    d_dt_theta = rps;
    
else 
    % ramp down
    d_dt_theta = max(-(rps/t_trans)*(t - t_down) + rps, 0);
    
end

% fprintf(1, '[shaftSpeed] time  = %g [s]\n', t);
% fprintf(1, '[shaftSpeed] theta = %g [rad]\n', theta);

end % shaftSpeed
