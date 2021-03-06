function VortexParticle(inputFile)
%VortexParticle.
%   VortexParticle(inputFile) executes the VortexParticle code for the case
%   described by the user parameters within the file inputFile. inputFile
%   is a string and is the short path to the filename (with the file extension included).

%% =======================================================================%
% Preamble
% ========================================================================%
close all;
fclose all;
clc;
profile off;
diary off;

% Add paths to the source code
if ~isdeployed
    % add folder to source code
    dir_src = pwd;
    addpath([dir_src filesep 'source'])
    addpath([dir_src filesep 'source' filesep 'PoissonFFT'])
    
    % add external dependencies
    dir_ext = [pwd filesep 'source' filesep 'external']; 
    addpath([dir_ext filesep 'bipolar_colormap'])
    addpath([dir_ext filesep 'consolidator'])
%     addpath([dir_ext filesep 'interparc'])
    addpath([dir_ext filesep 'plot3k'])
%     addpath([dir_ext filesep 'write_VTK'])
%     addpath([dir_ext filesep 'write_VTK' filesep 'vtk_writers'])
%     addpath([dir_ext filesep 'write_VTK' filesep 'vtk_writers' filesep 'VtkWriter-0.1'])
%     addpath([dir_ext filesep 'MOSAIC'])
end

% Set debug break points
% dbstop in induce_vel_gpu_v1.m at 1
% dbstop init_vortexRing.m 
dbstop in VortexParticle at 99
% dbstop in interp_M2P at 1
% dbstop in PoissonSolve3D at 63
% dbstop in vel_P2M_reg at 150
% dbstop if error

%% =======================================================================%
% Initialize part 1
% ========================================================================%
[SIM, ENV, MESH] = VortexParticle_Init(inputFile);

% if SIM.DEBUG_LVL > 0
%     % check the number of particles and ask for permission to continue
%     fprintf(1,'\n[DEBUG] The number of particles is %g, continue?\n', PART.nPart);
%     debugPrompt = input('[DEBUG] press enter to continue, or type ABORT to enter debug mode: ','s');
%     if strcmp(debugPrompt,'ABORT')
%         error('[DEBUG] ABORT')
%     end
% end

%% =======================================================================%
% Initialize part 2
% ========================================================================%                      
switch SIM.example
    
        case 'VIC'
        % @@@@@@@@@@@@@@@@@@@@@**^^""~~~"^@@^*@*@@**@@@@@@@@@
        % @@@@@@@@@@@@@*^^'"~   , - ' '; ,@@b. '  -e@@@@@@@@@
        % @@@@@@@@*^"~      . '     . ' ,@@@@(  e@*@@@@@@@@@@
        % @@@@@^~         .       .   ' @@@@@@, ~^@@@@@@@@@@@
        % @@@~ ,e**@@*e,  ,e**e, .    ' '@@@@@@e,  "*@@@@@'^@
        % @',e@@@@@@@@@@ e@@@@@@       ' '*@@@@@@    @@@'   0
        % @@@@@@@@@@@@@@@@@@@@@',e,     ;  ~^*^'    ;^~   ' 0
        % @@@@@@@@@@@@@@@^""^@@e@@@   .'           ,'   .'  @
        % @@@@@@@@@@@@@@'    '@@@@@ '         ,  ,e'  .    ;@
        % @@@@@@@@@@@@@' ,&&,  ^@*'     ,  .  i^"@e, ,e@e  @@
        % @@@@@@@@@@@@' ,@@@@,          ;  ,& !,,@@@e@@@@ e@@
        % @@@@@,~*@@*' ,@@@@@@e,   ',   e^~^@,   ~'@@@@@@,@@@
        % @@@@@@, ~" ,e@@@@@@@@@*e*@*  ,@e  @@""@e,,@@@@@@@@@
        % @@@@@@@@ee@@@@@@@@@@@@@@@" ,e@' ,e@' e@@@@@@@@@@@@@
        % @@@@@@@@@@@@@@@@@@@@@@@@" ,@" ,e@@e,,@@@@@@@@@@@@@@
        % @@@@@@@@@@@@@@@@@@@@@@@~ ,@@@,,0@@@@@@@@@@@@@@@@@@@
        % @@@@@@@@@@@@@@@@@@@@@@@@,,@@@@@@@@@@@@@@@@@@@@@@@@@    
        
            % The simulation algorithm thus proceeds as follows:
                % Initialize particles 
            % Do the time steps n = 0, . . . , T ? 1:
                % � interpolate ? from the particles onto the mesh
                % � solve ?v = ??�? for v on the mesh using a fast Poisson solver (e.g. multigrid or FFT
                % � interpolate v from the mesh back to the particles ? v
                % � compute vortex stretching ? � ?v on the mesh using, e.g., finite differences
                % � interpolate the vortex stretching result from the mesh to the particles as a
                % change of vorticity
                % � compute vorticity diffusion ? ?? (isotropic and homogeneous) on the particles
                % using PSE (or RW) (alternatively this can also be computed on the mesh),
                % add this vorticity change to the one from the stretching term, and update the
                % particle circulation. Also add the curl of the body force, ? � f
                % � move (advect) the particles with velocity v and update their positions, leaving the particle volumes unchanged since advection is incompressible.
                % � remesh if needed. Generate particles only at mesh nodes where |?| > 0.
            % end time step
        
           MESH           = init_mesh(SIM, MESH);               % init mesh
           wf             = init_field(SIM, MESH);              % init vorticity field


           [xp, wp, PART] = remesh_particles(SIM, MESH, wf);    % remesh the field onto new particles
%            wf             = interp_P2M(MESH, xp, wp);           % init vorticity field by interpolation from particles (P2M)
           wf2             = interp_P2M(MESH, xp, wp, 'timing');           % init vorticity field by interpolation from particles (P2M)
            dwf_x = abs( wf{1} - wf2{1} );
            dwf_y = abs( wf{2} - wf2{2} );
            dwf_z = abs( wf{3} - wf2{3} );
            dwf   = {dwf_x; dwf_y; dwf_z};
            plot_field(dwf, MESH, 'Vorticity Field')
            
           uf             = PoissonSolve3D(SIM, MESH, wf);      % solve Poisson eqn for velocity
           % compute vortex stretching on mesh (finite differences)
           wf_str = vortex_stretch(SIM, MESH, wf, uf);
           wp = wp + interp_M2P(MESH, xp, wf_str);              % interpolate vortex stretching from mesh to the particles as a change of vorticity
           % move (advect) the particles with the velocity and update thier positions (particle volumes unchanged since advection is incompressible)
           
           up         = interp_M2P(MESH, xp, uf);         % interpolate velocity field to particles (M2P)
           
           
           
           
           plot_field(wf, MESH, 'Vorticity Field')
           plot_field(uf, MESH, 'Velocity Field')
           plot_particles(xp, wp, MESH, 'Particle Vorticity')
           plot_particles(xp, up, MESH, 'Particle Velocity')
                   

        
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
        
   

        switch SIM.example
            case 'VortexRings'
                % particle positions, vorticity, and strengths.  
                % NOTE: is this initial field divergenge free (del dot omega)? should check this...
                [xp, wp, ap, PART.nPart] = init_vortexRing(RING, PART, MESH_SGS, ENV);

            case 'Turbine'
                % this type of simulation starts empty, and relies on "emitter" particles
                xp         = [];
                wp         = [];
                ap         = [];
                PART.nPart = 0;

            otherwise
        end

        %% =======================================================================%
        % compare the initial vorticity field to the remeshed particle approximation
        % ========================================================================%
        % recompute the field with interpP2M operation
        
        % recompute the field with full Biot-Savart law
        
        
        
         % initialize the mesh, extended mesh (for FFT), and vorticity field
        % [fieldk, fieldx_cen, fieldx_ext, fieldr_ext, NX, vort] = init_field(NXs, testcase, plot_InitField);
        [Mesh, wf] = init_field(SIM,MESH); % NOTE: Mesh is in NDGRID format, as opposed to the MESHGRID format
        
        
        
        
%         plot_field(wf, Mesh, 'Vorticity Field')

        
        % init particles with interpM2P()
%         [xp, ap] = init_particles(vort, Mesh, Sim)

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    case 'Turbine'       
        dt_emit     = 0.15; % emit particles every [seconds]
        
        theta     = 0; % initial value. shaft angular position.
        
        cycle     = 0;
        time      = 0;
        
        emit_t     = time:dt_emit:SIM.endtime;
        emit_cycle = 0;
        
        dt_out    = 1/SIM.fps_output;
        
        % pre-allocate all the particles
        xp     = [];
        ap     = [];
        xp_old = [];
        ap_old = [];
        
        xp_wake     = [];
        ap_wake     = [];
        xp_old_wake = [];
        ap_old_wake = [];
        wp_old_wake = [];
        up_old_wake = [];
        
        xp_inflow     = [];
        ap_inflow     = [];
        xp_old_inflow = [];
        ap_old_inflow = [];
        wp_old_inflow = [];
        up_old_inflow = [];
        
        while time < SIM.endtime-dt_out

            tspan  = [time, time+dt_out];
%             tspan = 0:0.1:SIM.endtime;
           
            theta_init  = theta;
            [T, Y]      = ode45(@shaftSpeed, tspan, theta_init, SIM.optionsODE, LL.ROT_SPD);
            theta_old   = Y(1); % solution at tspan(1)
            theta       = Y(2); % solution at tspan(end)
           
            azim_shaft_old  = mod(radtodeg(theta_old), 360);
            azim_shaft      = mod(radtodeg(theta),     360);
            azim_blades_old = mod(azim_shaft_old + 360/LL.NUM_BLADES .* (0:LL.NUM_BLADES-1)', 360);
            azim_blades     = mod(azim_shaft     + 360/LL.NUM_BLADES .* (0:LL.NUM_BLADES-1)', 360);
            d_dt_theta_old  = shaftSpeed(tspan(1), theta_old, LL.ROT_SPD);
            d_dt_theta      = shaftSpeed(tspan(2),     theta, LL.ROT_SPD);
            
            xi = cosspace(0, 1, LL.NUM_SEC, 'both'); % blade element spacing
            for n = 1:LL.NUM_BLADES
                % a structure array, containing transformations to each blade
                Coord(n)         = defineCoordSystems(LL, azim_blades(n));
                xp_emit_rotor{n} = interparc(xi, ...
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
            ap_rotor = [repmat(Coord(1).Tran.BG * [0; 0; 1], 1, LL.NUM_SEC), ...
                        repmat(Coord(2).Tran.BG * [0; 0; 1], 1, LL.NUM_SEC), ...
                        repmat(Coord(3).Tran.BG * [0; 0; 1], 1, LL.NUM_SEC)];
                        
            %% =======================================================================%
            % create particles subject to vorticity transport equation (VTE)
            % ========================================================================%            
            if tspan(2) > emit_t(emit_cycle+1)
                emit_cycle = emit_cycle + 1;
                
                % wake
                xp_emit_wake = xp_rotor;
                ap_emit_wake = [repmat(Coord(1).Tran.BG * [0; 0; 1], 1, LL.NUM_SEC), ...
                                repmat(Coord(2).Tran.BG * [0; 0; 1], 1, LL.NUM_SEC), ...
                                repmat(Coord(3).Tran.BG * [0; 0; 1], 1, LL.NUM_SEC)];

                % inflow turbulence
                [my, mz] = meshgrid(linspace(           -LL.ROTOR_DIA,             LL.ROTOR_DIA, 5), ...
                                    linspace(LL.HUB_HT - LL.ROTOR_DIA, LL.HUB_HT + LL.ROTOR_DIA, 5));
                xp_emit_inflow = [-LL.ROTOR_DIA/2.*ones(1,numel(my)); my(:)'; mz(:)'];
                ap_emit_inflow = rand(3, size(xp_emit_inflow, 2));
                
                % consolidate particles            
                xp        = [xp_old       , xp_emit_wake  , xp_emit_inflow];
                ap        = [ap_old       , ap_emit_wake  , ap_emit_inflow];
                xp_wake   = [xp_old_wake  , xp_emit_wake             ];
                ap_wake   = [ap_old_wake  , ap_emit_wake             ];
                xp_inflow = [xp_old_inflow, xp_emit_inflow           ];
                ap_inflow = [ap_old_inflow, ap_emit_inflow           ];
            else
                % consolidate particles
                xp        = xp_old;
                ap        = ap_old;
                xp_wake   = xp_old_wake;
                ap_wake   = ap_old_wake;
                xp_inflow = xp_old_inflow;
                ap_inflow = ap_old_inflow;
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
            [T, Y]  = ode113(@get_RHS_ode, tspan, x_init, SIM.optionsODE, PART, ENV, SIM);
            tmp     = reshape(Y(end,:)', 3*PART.nPart, 2);  % the solution matrix
            xp      = reshape( tmp(:,1),   PART.nPart, 3)';
            ap      = reshape( tmp(:,2),   PART.nPart, 3)';

            % instead use a push, pop, stack idea here
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
                MESH = updateMesh(MESH, LL.ROTOR_DIA/2, xp, PART);
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
            output_exampleTurbine(SIM, MESH, PART, ENV, Coord, LL, ...
                                  fig1, fig2, ...                          
                                  tspan, ...
                                  xp, ...
                                  xp_emit_rotor, ...
                                  azim_blades_old, ...
                                  azim_blades, ...
                                  d_dt_theta_old, ...
                                  d_dt_theta)
%             makeFigures_v2(MESH, LL, xp, ap, up, uf_x, uf_y, uf_z, wf_x, wf_y, wf_z, SIM.DEBUG_LVL)   
            
        end        
        
    case 'VortexRings'
        % =======================================================================%
        % update the mesh used for visualization
        % ========================================================================%
        if MESH.adaptive
            MESH = updateMesh(MESH, RING.Rmajor, xp, PART);
        end    
        
        % particle velocity (P2P)
        up = vel_P2P(xp, ap, PART, SIM, ENV);

        % RHS of the vorticity eqn (P2P)
        dap_dt_RHS = get_RHS(PART.nPart, xp, ap, PART.hp, ENV.kin_visc, SIM.runMode_RHS);

        % velocity field (P2M)
        [uf_x, uf_y, uf_z] = vel_P2M_reg(xp, ap, PART, SIM, MESH, ENV);

        % vorticity field (P2M)
        [wf_x, wf_y, wf_z] = vort_P2M(MESH, PART.nPart, xp, ap, PART.hp, SIM.runMode_P2M);

        % make figures for verification
        makeFigures(MESH, RING, xp, wp, up, uf_x, uf_y, uf_z, wf_x, wf_y, wf_z, Inf)
                
    otherwise
        
end


%% write output files of initial conditions
write_Point3D(PART.nPart, xp, wp, up, 0, SIM.outputDir)
save_VectorAndMagnitudeVTK_binary(uf_x,uf_y,uf_z,MESH.X,MESH.Y,MESH.Z, ...
                                  [SIM.outputDir filesep 'Velocity_Field_'  num2str(0,'%04.0f') '.vtk'],'vel_field')
save_VectorAndMagnitudeVTK_binary(wf_x,wf_y,wf_z,MESH.X,MESH.Y,MESH.Z, ...
                                  [SIM.outputDir filesep 'Vorticity_Field_' num2str(0,'%04.0f') '.vtk'],'vort_field')




%% find particles that are close enough to the boundary


%% create the image particles for this subset of particles


%% consolidate all particles into a single set


%% perform all the RHS summations (w/ new cutoff acceleration) and compute the velocity & vorticity fields


%% move the particles

%% check for any particles that leave the domain, or have vorticity below a cutoff value
% give remaining particle strength to nearby particles and remove from
% simulation (be careful about unneccesary resizing of array each timestep)




%% =======================================================================%
% Solver: Begin time stepping the evolution equations.
%         At this point, need to be given initial: 
%         mesh, xp, ap (wp), up, uf, wf
% ========================================================================%
if SIM.DEBUG_LVL > 0
    % ask permission to start main solver
    fprintf(1,'\n[DEBUG] Does everything look good...should we start the main solver?\n',PART.nPart);
    debugPrompt = input('[DEBUG] press enter to continue, or type ABORT to do that thing...','s');
    if strcmp(debugPrompt,'ABORT')
        error('[DEBUG] ABORT')
    end
end
fprintf(['\nInitialization complete. Executing for case ' SIM.caseName '.\n'])

% set the profiling settings
if SIM.DEBUG_LVL > 9000
    profile on
end


%% new time stepping loop
cycle  = 0;
time   = 0;
dt_out = 1/SIM.fps_output;   
while time < SIM.endtime-dt_out

    tspan   = [time, time+dt_out];

    x_init  = [reshape(xp',3*PART.nPart,1); 
               reshape(ap',3*PART.nPart,1)];

%     [T, Y]  = ode45(@get_RHS_ode, tspan, x_init, SIM.optionsODE, PART, ENV, SIM);
    [T, Y]  = ode113(@get_RHS_ode, tspan, x_init, SIM.optionsODE, PART, ENV, SIM);
    
    tmp = reshape(Y(end,:)', 3*PART.nPart, 2);  % the solution matrix
    xp  = reshape( tmp(:,1),   PART.nPart, 3)';
    ap  = reshape( tmp(:,2),   PART.nPart, 3)';

    %% =======================================================================%
    % update particle velocities: velocity at particle locations
    % ========================================================================%
    up = vel_P2P(xp, ap, PART, SIM, ENV);    %NOTE: this is already computed in the get_RHS_ode function!

    if true                             %NOTE: add a better condition for this
        %% =======================================================================%
        % update the mesh used for visualization
        % ========================================================================%
        if MESH.adaptive
            MESH = updateMesh(MESH, RING, xp, PART);
        end

        %% =======================================================================%
        % update velocity field: velocity at mesh point locations
        % ========================================================================%
        [uf_x, uf_y, uf_z] = vel_P2M_reg(MESH, PART.nPart, xp, ap, PART.hp, SIM.runMode_P2M, ENV);

        %% =======================================================================%
        % update vorticity field: vorticity at mesh point locations
        % ========================================================================%
        [wf_x, wf_y, wf_z] = vort_P2M(MESH, PART.nPart, xp, ap, PART.hp, SIM.runMode_P2M);

        %% =======================================================================%
        % write the output files
        % ========================================================================%
        if SIM.writeParticles
            write_Point3D(PART.nPart, xp, ap, up, cycle, SIM.outputDir)
        end
        if SIM.writeVelocityField
            save_VectorAndMagnitudeVTK_binary(uf_x, uf_y, uf_z, MESH.X, MESH.Y, MESH.Z, [SIM.outputDir filesep 'Velocity_Field_'  num2str(cycle,'%04.0f') '.vtk'],'vel_field')
        end
        if SIM.writeVorticityField
            save_VectorAndMagnitudeVTK_binary(wf_x, wf_y, wf_z, MESH.X, MESH.Y, MESH.Z, [SIM.outputDir filesep 'Vorticity_Field_' num2str(cycle,'%04.0f') '.vtk'],'vort_field')
        end
        makeFigures(MESH, RING, xp, ap, up, uf_x, uf_y, uf_z, wf_x, wf_y, wf_z, SIM.DEBUG_LVL)
        
        %% =======================================================================%
        % This is a good place to check other diagnostics
        % ========================================================================%
        % NOTE: a better condition to terminate would be when the particle cores no longer overlap! 
        if assignParticleMask(MESH, xp)
            fprintf(1, '[ERROR] particles are outside the mesh domain! \n');
            die
        end
        if any( [MESH.dx, MESH.dy, MESH.dz] > MESH.max_dh )
            MESH
            fprintf(1, '[ERROR] the adaptive mesh has grown too large! \n');
            die
        end
        
    end 
    
    % update the time span for next integration cycle
    cycle = cycle + 1;
    time  = tspan(2);

end

%% =======================================================================%
% final saving of data and cleanup
% ========================================================================%
save([SIM.outputDir filesep SIM.caseName '_workspace.mat'])
fprintf(['All data in workspace saved to MAT-file: ' SIM.caseName '_workspace.mat.\n'])
fprintf([SIM.version ' terminated normally for case ' SIM.caseName '.\n'])
diary OFF

if SIM.DEBUG_LVL > 9000
    % profile the code
    profile viewer
end

% %% tear down parallel computing stuff
% matlabpool close force local


end % function VortexParticle

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