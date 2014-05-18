% Needs from Computation/Theory 
% * Dynamic mesh adaptation to help resolve/preserve key flow structures/processes and bound computational workload 
% * Adaptative solver algorithms that span the DNS/LES/DES/RANS hierarchy to improve physical fidelity and manage computational workload 
% * Tools to interrogate extremely large computational and experimental data sets, in a way that exploits synergies between predicted and measured data 

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
% @@ here be dragons @@@@@~ ,@@@,,0@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@,,@@@@@@@@@@@@@@@@@@@@@@@@@

%% =======================================================================%
% Simulation Parameters
% ========================================================================%
SIM.outputDir           = 'C:\Users\Danny\Desktop\simulation_output\METS-2014\test_VIC-actuator-lines';
% SIM.outputDir           = '/home/danny/workspace/simulation_output/vic-turbine-1';
SIM.DEBUG_LVL           = 8999;   	% setting a debug level > 0 shows additional output.  If you go over 9000 the profiler is enabled.
SIM.writeParticles      = true;
SIM.writeVelocityField  = true;
SIM.writeVorticityField = true;
SIM.writePressureField  = true;

%==========================================================================
% Input options 
%==========================================================================
%% settings for the FFT Poisson solver
% kernel: Regularisation order of integration kernel 
%   0  = non-regularised
%   2  = 2nd order
%   4  = 4th order
%   6  = 6th order
%   8  = 8th order
%   10 = 10th order
% solve_vel: Type of Poisson equation
%   0 = solve for stream function by G kernel
%   1 = solve for velocity by K kernels, 
%   2 = solve for velocity by G kernel + spectral differentiating (only works for periodic BC, or also freespace BC?)
% alpha: Smoothing radius relative to mesh size: epsilon = alpha*dx (default 2) (and where is this used?)
SIM.kernel    = 4;
SIM.solve_vel = 1; 
SIM.alpha     = 2; 
SIM.domainbc  = 0; % freespace boundary conditions (=1 for periodic - not fully tested yet)

%% algorithm selection for P2P, P2M, and RHS evaluations -- these are mostly for the old PSE routines
SIM.numProc     = 4;        % number of CPU cores to use for parallel algorithms
SIM.runMode_P2P = 'GPU-v2'; % choose: 'CPU-v1', 'CPU-v2', 'CPU-v3', 'GPU-v1', 'GPU-v2'
SIM.runMode_RHS = 'GPU-v2'; % choose: 'CPU-v1', 'CPU-v2', 'CPU-v3', 'GPU-v1', 'GPU-v2'
SIM.runMode_P2M = 'GPU-v2'; % choose: 'CPU-v1', 'CPU-v2',           'GPU-v1', 'GPU-v2'

%% set time-stepping and output frequency
SIM.endtime    = 5;
SIM.fps_output = 10;
SIM.TSTEPPING  = 'fixed';    % choose 'fixed' or 'variable' time stepping algorithms
SIM.dt         = 0.1;   % time step (only used for fixed step integrators)
SIM.optionsODE = odeset('AbsTol',           1e-4, ...
                        'RelTol',           1e-4, ...
                        'MaxStep',          1/SIM.fps_output, ...
                        'InitialStep',      [], ...
                        'NormControl',      'on', ...
                        'NonNegative',      [], ...
                        'Refine',           1, ...
                        'Stats',            'off', ...
                        'Mass',             [], ...
                        'MStateDependence', 'none', ...
                        'Events',           [], ...
                        'OutputFcn',        [], ...
                        'OutputSel',        []);    %NOTE: these options only used by Matlab ODE solvers, not 3rd party solvers

%% =======================================================================%
% Environmental & Fluid Properties
% ========================================================================%
% ENV.kin_visc = 1.46e-5;         % fluid kinematic viscosity (m^2/s) AIR
ENV.kin_visc = 1.05e-6;      	% fluid kinematic viscosity (m^2/s) WATER
% ENV.velFree  = [0; 0; 0];   	% Free stream velocity (a 3x1 array) [m/s]
ENV.velFree  = [1; 0; 0];   	% Free stream velocity (a 3x1 array) [m/s]

%% =======================================================================%
% Particle & Mesh Parameters
% ========================================================================%
SIM.dim        = 3;       % spatial dimensions (2D and 3D supported)
SIM.mbc        = 2;       % ghost layer, size of finite diff stencil beyond boundaries
SIM.h_cutoff   = 2;       % (same as alpha?) used to define support of a particles in terms of the mesh spacing (hp = h_cutoff * mesh.dx).  stay within range hp/dx > 1
SIM.tol_remesh = 1e-6;    % set vorticity to zero when the field is less than cutoff percent of the field maximum
SIM.pad        = 5;       % minimum distance between mesh boundaries and particles w.r.t. particle support (MESH.pad = SIM.pad * PART.hp;)

% parent mesh
MESH.tag      = 'parent mesh';
MESH.type     = 'colocated';    % staggered meshes are not yet supported
MESH.xmin     = 2*[-1, -1, -1];
MESH.xmax     = 2*[ 1,  1,  1];
MESH.NX       = [24, 24, 24];  % use powers of 2^nx
MESH.adaptive = false;

%% =======================================================================%
% Example Specific Parameters
% ========================================================================%
CTRL.testcase = 4;
switch CTRL.testcase
    case 1
        %% Bump function: SPHERICAL SCALAR FIELD
        CTRL.param.c = 20; % Function constant
        CTRL.param.R = 1;  % Function constant
        
    case 2
        %% Bump function: VORTEX RING (in the xy-plane)
        CTRL.param.c = 10;   % Function constant
        CTRL.param.R = 1;    % Function constant
        
    case 3
       %% =======================================================================%
        %   _   _            _             ______ _                 
        %  | | | |          | |            | ___ (_)                
        %  | | | | ___  _ __| |_ _____  __ | |_/ /_ _ __   __ _ ___ 
        %  | | | |/ _ \| '__| __/ _ \ \/ / |    /| | '_ \ / _` / __|
        %  \ \_/ / (_) | |  | ||  __/>  <  | |\ \| | | | | (_| \__ \
        %   \___/ \___/|_|   \__\___/_/\_\ \_| \_|_|_| |_|\__, |___/
        %                                                  __/ |    
        %                                                 |___/ 
        % ========================================================================%
        % need to specify for each vortex ring:
        CTRL.Re       = [10000];                         % Reynolds number of the vortex ring, defined as Re = gamma / kin_visc [ring1, ring2, ...]
        CTRL.Rmajor   = [1];                            % major radius of vortex ring [ring1, ring2, ...]
        CTRL.Rminor   = 0.2 .* CTRL.Rmajor;             % minor radius of vortex ring [ring1, ring2, ...]
        CTRL.center_x = [0];                            % x-coordinate of ring center [ring1, ring2, ...]
        CTRL.center_y = [0];                            % y-coordinate of ring center [ring1, ring2, ...]
        CTRL.center_z = [0];                            % z-coordinate of ring center [ring1, ring2, ...]
        CTRL.sign     = [1];
        % not yet implemented:
        % CTRL.axis_x   = [0, 0];                      	% x-component of unit vector defining the axis which the vortex ring is aligned [ring1, ring2, ...]
        % CTRL.axis_y   = [0, 0];                     	% y-component of unit vector defining the axis which the vortex ring is aligned [ring1, ring2, ...]
        % CTRL.axis_z   = [1, -1];                    	% z-component of unit vector defining the axis which the vortex ring is aligned [ring1, ring2, ...]
        % CTRL.azimAmp  = [0, 0];                       % amplitude of azimuthal purturbation to the circulation strength [ring1, ring2, ...]
        % CTRL.azimFreq = [0, 0];                       % frequency of azimuthal purturbation to the circulation strength [ring1, ring2, ...]
        % CTRL.radAmp   = [0, 0];                       % amplitude of radial purturbation to the circulation strength [ring1, ring2, ...]
        % CTRL.radFreq  = [0, 0];                       % frequency of radial purturbation to the circulation strength [ring1, ring2, ...]
        
    case 4
       %% =======================================================================%
        %   _____          _     _            
        %  |_   _|        | |   (_)           
        %    | |_   _ _ __| |__  _ _ __   ___ 
        %    | | | | | '__| '_ \| | '_ \ / _ \
        %    | | |_| | |  | |_) | | | | |  __/
        %    \_/\__,_|_|  |_.__/|_|_| |_|\___|
        % ========================================================================%
        % lifting line simulation
        CTRL.NUM_BLADES  = 3;            % Number of blades
        CTRL.ROTOR_DIA   = 20;           % Rotor diameter [m]
        CTRL.HUB_DIA     = 2;            % Hub diameter [m]
        CTRL.HUB_HT      = 20;           % Hub height [m]
        CTRL.NUM_SEC     = 30;           % Number of blade cross sections
        CTRL.ROT_SPD     = 11.5;         % Rotor rotational speed [rpm]
        CTRL.HH_SPD      = 1.0;          % Free stream flow speed [m/s] at hub height (hub height = z distance to center of the energy extraction area)
        CTRL.YAW         = 0;
        CTRL.SHAFT_TILT  = 0;
        CTRL.PRE_CONE    = 0;
        CTRL.TEETER      = 0;
        CTRL.BLD_PITCH   = 0;
        CTRL.INFLOW_TURB = false;
        % for now, some other parameter are hard coded in for rotor geometry and controls (DOE Ref. Model 1 - Tidal Turbine - But could read input files from either WT_Perf/FAST/CACTUS ??? can we choose some standard like in OpenMDAO ???)

        
end



