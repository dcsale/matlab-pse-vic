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
% dbstop in case_vortexRings at 19
% dbstop in PoissonSolve3D at 38
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
switch CTRL.testcase        
    case {2,3}
        case_vortexRings(CTRL, SIM, MESH, ENV)       
    otherwise
        error('[ERROR] testcase not recognized')       
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

