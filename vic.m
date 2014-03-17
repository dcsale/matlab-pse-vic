function vic(inputFile)

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
    addpath([dir_ext filesep 'interparc'])
    addpath([dir_ext filesep 'plot3k'])
    addpath([dir_ext filesep 'ODE_Solvers'])
%     addpath([dir_ext filesep 'write_VTK'])
    addpath([dir_ext filesep 'write_VTK' filesep 'vtk_writers'])
%     addpath([dir_ext filesep 'write_VTK' filesep 'vtk_writers' filesep 'VtkWriter-0.1'])
%     addpath([dir_ext filesep 'MOSAIC'])
end

% Set debug break points
% dbstop in vic at 110
% dbstop in case_vortexRings at 1
dbstop if error

%% =======================================================================%
% Initialize part 1
% ========================================================================%
[SIM, ENV, MESH, CTRL] = VortexParticle_Init(inputFile);

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
    case 'VortexRings'
        case_vortexRings(CTRL, SIM, MESH, ENV)       
    case 'Turbine'          
        case_turbine(CTRL, SIM, MESH, ENV)
        fprintf(1, 'done.\n');
        return
    otherwise
        error('SIM.example not recognized')       
end


%% write output files of initial conditions
write_Point3D(PART.nPart, xp, wp, up, 0, SIM.outputDir)
save_VectorAndMagnitudeVTK_binary(uf_x,uf_y,uf_z,MESH.X,MESH.Y,MESH.Z, ...
                                  [SIM.outputDir filesep 'Velocity_Field_'  num2str(0,'%04.0f') '.vtk'],'vel_field')
save_VectorAndMagnitudeVTK_binary(wf_x,wf_y,wf_z,MESH.X,MESH.Y,MESH.Z, ...
                                  [SIM.outputDir filesep 'Vorticity_Field_' num2str(0,'%04.0f') '.vtk'],'vort_field')





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

% %% tear down parallel computing stuff
% matlabpool close force local

if SIM.DEBUG_LVL > 9000
    % profile the code
    profile viewer
end


end % function VortexParticle

