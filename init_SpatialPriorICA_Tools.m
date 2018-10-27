BaseDir = fileparts(mfilename('fullpath'));
addpath([BaseDir,filesep,'ClusterAnalysisTools']);
disp(['Cluster Analysis tools have been added. ("',BaseDir,filesep,'ClusterAnalysisTools','")']);
addpath([BaseDir,filesep,'ParcellationTools']);
disp(['Parcellation tools have been added. ("',BaseDir,filesep,'ParcellationTools','")']);
addpath([BaseDir,filesep,'DisplayROI_ICs']);
disp(['Display tools for CompEncode Overlay generation have been added. ("',BaseDir,filesep,'DisplayROI_ICs','")']);
addpath([BaseDir,filesep,'ClusterTools']);
disp(['Clustering tools (taken from AtlasROITools V2.0) have been added. ("',BaseDir,filesep,'ClusterTools','")']);
init_ClusterTools;