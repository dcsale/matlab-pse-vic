% Needs from Computation/Theory 
% ? Dynamic mesh adaptation to help resolve/preserve key flow structures/processes and bound 
% computational workload 
% ? Adaptative solver algorithms that span the DNS/LES/DES/RANS hierarchy to improve 
% physical fidelity and manage computational workload 
% ? Algorithms to interrogate extremely large computational and experimental data sets, in a way 
% that exploits synergies between predicted and measured data. 

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

%% =======================================================================%
% Simulation Parameters
% ========================================================================%
% SIM.outputDir           = 'C:\Users\Danny\Desktop\simulation_output\METS-2014\VortexInCell - test - 1';
SIM.outputDir           = '/home/danny/workspace/simulation_output/VortexInCell-test-1';
SIM.example             = 'VIC';                % current options are: 'VortexRings', 'Turbine', "VIC'
SIM.DEBUG_LVL           = 8999;                 % setting a debug level > 0 shows additional output.  If you go over 9000 the profiler is enabled.
SIM.writeParticles      = true;
SIM.writeVelocityField  = false;
SIM.writeVorticityField = false;

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
SIM.kernel    = 4;

% solve_vel: Type of Poisson equation
%   0 = solve for stream function by G kernel
%   1 = solve for velocity by K kernels, 
%   2 = solve for velocity by G kernel + spectral differentiating
SIM.solve_vel = 1; 

% alpha: Smoothing radius relative to mesh size: epsilon = alpha*dx (default 2)
SIM.alpha     = 2; 

% NXs: Number of mesh cells (use vector for convergence studies)
% Sim.NXs       = 64*2.^(0:2);
% SIM.NXs       = 8;

% Testcases
% 1 = Bump function: SPHERICAL SCALAR FIELD (only for solve_vel = 0)
% 2 = Bump function: VORTEX RING (xy-plane)
SIM.testcase = 2;

%% algorithm selection for P2P, P2M, and RHS evaluations
SIM.numProc     = 4;        % number of CPU cores to use for parallel algorithms
SIM.runMode_P2P = 'GPU-v2'; % choose: 'CPU-v1', 'CPU-v2', 'CPU-v3', 'GPU-v1', 'GPU-v2'
SIM.runMode_RHS = 'GPU-v2'; % choose: 'CPU-v1', 'CPU-v2', 'CPU-v3', 'GPU-v1', 'GPU-v2'
SIM.runMode_P2M = 'GPU-v2'; % choose: 'CPU-v1', 'CPU-v2',           'GPU-v1', 'GPU-v2'
%% set time-stepping and output frequency
SIM.endtime    = 5;
SIM.fps_output = 15;
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
                        'OutputSel',        []);

%% =======================================================================%
% Environmental & Fluid Properties
% ========================================================================%
% ENV.kin_visc = 1.46e-5;       % fluid kinematic viscosity (m^2/s) AIR
ENV.kin_visc = 1.05e-6;       % fluid kinematic viscosity (m^2/s) WATER
ENV.velFree  = [1; 0; 0];           % Free stream velocity (a 3x1 array) [m/s]

%% =======================================================================%
% Particle & Mesh Parameters
% ========================================================================%
% vortex particles
PART.h_cutoff     = 1.5;        % used to define support of a single vortex particle in terms of the mesh spacing (hp = h_cutoff * mesh.dx).  stay within range hp/dx > 1
% PART.hp           = 0.2;        % a smoothing radius (i.e., a cutoff length or core size)

% parent mesh
MESH.type         = '3D cartesian: rectilinear';
MESH.adaptive     = true;
MESH.NX           = 2*[9, 9, 9];
% "sub-grid scale" mesh - this mesh is only used for initialization and remeshing of vortex particles
% MESH_SGS.type     = '3D cartesian: rectilinear';
% MESH_SGS.adaptive = true;
% MESH_SGS.NX       = [128; 128; 128];

%% =======================================================================%
% Example Specific Parameters
% ========================================================================%
if SIM.testcase == 1 
    % Bump function: SPHERICAL SCALAR FIELD
    SIM.param.c = 20; % Function constant
    SIM.param.R = 1;  % Function constant

elseif SIM.testcase == 2 
    % Bump function: VORTEX RING (in the xy-plane)
    SIM.param.c = 10;  % Function constants
    SIM.param.R = 0.5; % Function constants

end
%% =======================================================================%
%   _____          _     _            
%  |_   _|        | |   (_)           
%    | |_   _ _ __| |__  _ _ __   ___ 
%    | | | | | '__| '_ \| | '_ \ / _ \
%    | | |_| | |  | |_) | | | | |  __/
%    \_/\__,_|_|  |_.__/|_|_| |_|\___|
% ========================================================================%
%% lifting line simulation
LL.NUM_BLADES = 3;            % Number of blades
LL.ROTOR_DIA  = 20;           % Rotor diameter [m]
LL.HUB_DIA    = 2;            % Hub diameter [m]
LL.HUB_HT     = 20;           % Hub height [m]
LL.NUM_SEC    = 10;           % Number of blade cross sections
LL.ROT_SPD    = 11.5;         % Rotor rotational speed [rpm]
LL.HH_SPD     = 3.5;          % Free stream flow speed [m/s] at hub height (hub height = z distance to center of the energy extraction area)
LL.YAW        = 30;
LL.SHAFT_TILT = 0;
LL.PRE_CONE   = 0;
LL.TEETER     = 0;
LL.BLD_PITCH  = 0;
% for now, some other parameter are hard coded in (DOE Ref. Model 1 - Tidal Turbine - But could read input files from WT_Perf/FAST)

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
RING.Re       = [3000];                  	% Reynolds number of the vortex ring, defined as Re = gamma / kin_visc [ring1, ring2, ...]
RING.Rmajor   = [1];                         % major radius of vortex ring [ring1, ring2, ...]
RING.Rminor   = 0.05 .* RING.Rmajor;           	% minor radius of vortex ring [ring1, ring2, ...]
RING.center_x = [0];                         % x-coordinate of ring center [ring1, ring2, ...]
RING.center_y = [0];                         % y-coordinate of ring center [ring1, ring2, ...]
RING.center_z = [1];                         % z-coordinate of ring center [ring1, ring2, ...]
RING.sign     = [1];
% not yet implemented:
% RING.axis_x   = [0, 0];                      	% x-component of unit vector defining the axis which the vortex ring is aligned [ring1, ring2, ...]
% RING.axis_y   = [0, 0];                     	% y-component of unit vector defining the axis which the vortex ring is aligned [ring1, ring2, ...]
% RING.axis_z   = [1, -1];                    	% z-component of unit vector defining the axis which the vortex ring is aligned [ring1, ring2, ...]
% RING.azimAmp  = [0, 0];                       % amplitude of azimuthal purturbation to the circulation strength [ring1, ring2, ...]
% RING.azimFreq = [0, 0];                       % frequency of azimuthal purturbation to the circulation strength [ring1, ring2, ...]
% RING.radAmp   = [0, 0];                       % amplitude of radial purturbation to the circulation strength [ring1, ring2, ...]
% RING.radFreq  = [0, 0];                       % frequency of radial purturbation to the circulation strength [ring1, ring2, ...]
