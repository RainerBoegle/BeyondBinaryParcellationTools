function [CompEncodeOutDir,CompEncodeOutPaths,MELODICdir,MMthresh] = CreateCompEncodeFromMELODICresults(MELODICdir,MMthresh,UseReThresh,Selection)
% This function creates the component encode from MELODIC results, e.g. saved in "groupmelodic.ica"-directory.
%
% Created will be a SimpleEncode with all voxels that show overlaps set to "-1" and pure/singular areas get the component index,
% a Base2Encode with unique areas set to powers of 2 and overlaps indicated by sums of powers of 2. (i.e. clear association via the sum, e.g. 3 is 2+1; 5 is 4+1; 7 is 4+2+1...),
% a Log2Base2Encode, i.e. just the log2 of the Base2Encode (>0) PLUS 1 such that a display with (rough) distinction of overlaps can be made.
% a SumMap, i.e. the summation of all significant voxels, indicating per voxel how many significant voxels over all compontent are there.
% a WinnerTakeAllEncode, i.e. simply checking all components for the highest z-stats value and assigning this component at the voxel. (Simplified version of the selection in the MultiLevelICA paper by Kim et al HBM 2012.)
% a mask ([0,1]) indicating all OverlapsOnly.
% a UniqueEncode that indicates the pure7singular areas by component index, i.e the simple encode without all areas that are "-1".
% 
%
%
%Usage:
%      [CompEncodeOutDir,CompEncodeOutPaths,MELODICdir,MMthresh] = CreateCompEncodeFromMELODICresults(MELODICdir,MMthresh,UseReThresh,Selection);
%      [CompEncodeOutDir,CompEncodeOutPaths,MELODICdir,MMthresh] = CreateCompEncodeFromMELODICresults(); %manual input/selection
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (25.08.2015): initial implementation

%% check inputs
%CATCH_MELODICdir
try
    if(~exist(MELODICdir,'dir'))
        error(['MELODICdir: "',MELODICdir,'" not found!']);
    else
        if(~exist([MELODICdir,filesep,'melodic_IC.nii'],'file'))
            error(['"melodic_IC.nii" not found in MELODICdir: "',MELODICdir,'"!']);
        elseif(~exist([MELODICdir,filesep,'melodic_mix'],'file'))
            error(['"melodic_mix" not found in MELODICdir: "',MELODICdir,'"!']);
        elseif(~exist([MELODICdir,filesep,'mask.nii'],'file'))
            disp(['"mask.nii" not found in MELODICdir: "',MELODICdir,'"!']);
            disp( 'No problem, will do mixture model thresholding without mask, ie use implicit mask.');
            UseMask = 0;
        end
    end
catch CATCH_MELODICdir
    MELODICdir = spm_select(1,'dir','Select MELODICdir, ie the directory containing the results of melodic for fMRI data...');
    disp_catch(CATCH_MELODICdir,[mfilename,'>CATCH_MELODICdir'],'CATCH_MELODICdir'); %assignin('base','CATCH_MELODICdir',CATCH_MELODICdir);
end
%CATCH_MMthresh  
try
    if(isempty(MMthresh))
        disp('MMthresh is empty! Setting it to 0.5 (default).');
        MMthresh = 0.5;
    else
        if(~((MMthresh>0)&&(MMthresh<1))) %MMthresh has to be between 0&1!!!
            disp('MMthresh has to be between 0&1!!! Setting it to 0.5 (default).');
            MMthresh = 0.5;
        end
    end
catch CATCH_MMthresh
    ans_MMthresh = inputdlg({'MixtureModel threshold p= '},'Threshold?',1,{'0.5'});
    disp_catch(CATCH_MMthresh,[mfilename,'>CATCH_MMthresh'],'CATCH_MMthresh'); %assignin('base','CATCH_MMthresh',CATCH_MMthresh);
    MMthresh = eval(ans_MMthresh{1});
end
%CATCH_UseReThresh
try
    if(isempty(UseReThresh))
        disp('"UseReThresh" not set! Will default to NOT using ReThresholding of stats-maps produced by MixtureModel-thresholding.');
        UseReThresh = 0;
    else
        if(UseReThresh<=0)
            UseReThresh = 0;
            disp('Will NOT use ReThresholding of stats-maps produced by MixtureModel-thresholding!');
        elseif(UseReThresh>0)
            UseReThresh = 1;
            disp('Will USE ReThresholding of stats-maps produced by MixtureModel-thresholding!');
        end
    end
catch CATCH_UseReThresh
    disp_catch(CATCH_UseReThresh,[mfilename,'>CATCH_UseReThresh'],'CATCH_UseReThresh'); %assignin('base','CATCH_Selection',CATCH_Selection);
    disp('"UseReThresh" not set! Will default to NOT using ReThresholding of stats-maps produced by MixtureModel-thresholding.');
    UseReThresh = 0;
end
%CATCH_Selection
try
    if(~isempty(Selection))
        if(isinf(Selection))
            disp('Using all ICs...');
        end
    end
catch CATCH_Selection
    Selection = []; %user select later
    disp_catch(CATCH_Selection,[mfilename,'>CATCH_Selection'],'CATCH_Selection'); %assignin('base','CATCH_Selection',CATCH_Selection);
end

%% create mixture model outputs
MixModelOutDir = [MELODICdir,filesep,'MixModelThreshResults_',datestr(now,'yyyymmmdd_HHMM')];
[MixModelOutDir,MMthresh,MELODICdir] = CreateMixModelThreshMaps(MELODICdir,MMthresh,MixModelOutDir);

%% ReThresh?
if(UseReThresh)
    answer_ReThresh = inputdlg({'ReThresh-PercentileThreshold= ';'ClusterSizeThreshold= '},'ReThresh Inputs',1,{'0.5';'13'});
    ThreshPercentile= eval(answer_ReThresh{1});
    ClusterSizeThres= eval(answer_ReThresh{2});
    
    files = cellstr(spm_select('List','^thresh_zstat.*',[MixModelOutDir,filesep,'stats']));
    if(isempty(files))
        error(['Could not find "^thresh_zstat.*" NIFTI-files in Directory "',MixModelOutDir,filesep,'stats".']);
    else
        for IndFile = 1:length(files)
            ReThreshICmaps(files{IndFile},ThreshPercentile,ClusterSizeThres);
        end
    end
end

%% create CompEncode
MixModelStatsDir = [MixModelOutDir,filesep,'stats'];
[CompEncodeOutDir,CompEncodeOutPaths] = CreateCompEncodeFromThreshMaps(MixModelStatsDir,UseReThresh,Selection);
save([CompEncodeOutDir,filesep,'MMthreshSettings.mat'],'MMthresh','MixModelOutDir','MELODICdir');

%% done.
disp('Done creating CompEncode from MELODIC results.');

end