function [SIM, ENV, MESH] = VortexParticle_Init(inputFile)

%% set program version
SIM.version = 'vic-alpha';

%% read/evaluate the input file variables into workspace
[~, SIM.caseName, SIM.inpExt] = fileparts(inputFile);
run(inputFile)
       
%% Define the files in the solution path
SIM.rootDir             = pwd;
SIM.inpFile             = [SIM.rootDir filesep SIM.caseName SIM.inpExt];
SIM.logfile             = [SIM.outputDir filesep SIM.caseName '_Log.txt'];
% SIM.xlsfile             = [SIM.outputDir filesep SIM.caseName '.xls'];
% SIM.airfoilInputData    = [SIM.rootDir   filesep 'Airfoil_Data'];
% SIM.optimizationOptions = [SIM.rootDir   filesep 'Optimization_Options'];
   
%% setup output directory 
if exist(SIM.outputDir,'dir') == 7
    debugPrompt = input('[WARNING] output directory exists and may contain important data! Should we clear the data? [yes / no]: ','s');
    if strcmp(debugPrompt,'yes')
        rmdir([SIM.outputDir],'s');
        mkdir(SIM.outputDir);
    else
        error('[ABORT]: you should modify the variable SIM.outputDir to avoid this error message again.');
    end  
else
    mkdir(SIM.outputDir);
end

%% echo user parameters to output directory
copyfile([SIM.inpFile],[SIM.outputDir filesep SIM.caseName '_echo' SIM.inpExt])

%% Setup parallel computing stuff
% if strcmp(SIM.runMode_P2P,'CPU-v1') || strcmp(SIM.runMode_P2M,'CPU-v1')
if SIM.numProc <= 1  
    matlabpool close force local
else
    poolSize = matlabpool('size');  % check to see if a pool is already open
    if poolSize == 0 || poolSize < SIM.numProc || poolSize ~= SIM.numProc
        matlabpool close force local
        eval(['matlabpool open ' num2str(SIM.numProc)])
    end
end

%% Start a log of all text and command line output. Read this file if
% something goes wrong and follow the error handling messages to debug.
diary(SIM.logfile);

end % function VortexParticle_Init
