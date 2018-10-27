function [MapExtractStruct,xSPM] = GetParamsFromSPMmat(SPMmatDir,SPMthresStruct)
% This function can extract the statistic values & coordinates from a contrast statistics evaluation in SPM-results.
% This can be based on a SPM.mat file, contrast selection and statistics parameters in the SPMthresStruct,
% OR it can be entered manually.
%
%INPUTS:
%       SPMmatDir                   <--     Directory containing SPM.mat file. Should be a string (If empty, then manual select via spm_select.)
%       SPMthresStruct.
%                     .Contrast     <--     For selecting the contrast. This should either be a string or a number. If string then it should fit the name of the desired contrast; If number then it indicates the contrast via this number.';
%                                           NB: If it is left empty (i.e. [] or missing), then you will be asked to select it manually
%                     .MultTestCor  <--     The multiple correction method to be used. Should be a string, either ''FWE'' or ''none''. (single quotations of course!)'; %"FDR" might work, but I doubt it, -not checked.
%                                           NB: If it is left empty (i.e. [] or missing), then you will be asked to select it manually
%                     .p            <--     The threshold p-value. Should be a number, e.g. 0.05 or 0.001 or such.';
%                                           NB:
%                                              If it is left empty (i.e. [] or missing), then SPMthresStruct.MultTest_Cor will be checked
%                                              if 'FWE'     then SPMthresStruct.p= 0.05;
%                                              if 'none'    then SPMthresStruct.p= 0.001;
%                                              if [](empty) then SPMthresStruct.p= []; (empty) and you will be asked to select it manually
%                     .k            <--     The minimum cluster size in voxels. Should be a number (integer>=0),
%                                           e.g. 0 or 6 or 27 or 42 or 81 or such.';                 
%                                           NB: If it is left empty (i.e. [] or missing), then you will be asked to select it manually
%OUTPUTS:
%       MapExtractStruct.
%                       .Thresholds   <--   Positive & Negative threshold for Map.
%                                           NB: NOT A P-VALUE! This will be a t- or F-value (i.e. xSPM.u) given the statistics settings in the SPM.mat and the contrast.
%
%                       .Negative.    <--   Negative threshold data (NB: HERE ALWAYS EMPTY, i.e. not used.)
%                       .Positive.    <--   Positive threshold data
%                                .Coords_mm   <<   The mm-coordinates of all significant voxels, i.e. "above" threshold.
%                                .StatsVals   <<   The statistics values of all significant voxels, i.e. "above" threshold.
%                                .Coords_vox  <<   The voxel-coordinates of all significant voxels, i.e. "above" threshold.
%                                .ConNr       <<   The number of the contrast
%
%                       .V_map.       <--   SPM-vol struct of the input map (can be used for later output of NIFTI).
%                                           NB: this will point to the unthresholded SPM NIFTI with the t- or F-values.
%
%Usage:
%       [MapExtractStruct] = GetParamsFromSPMmat(SPMmatDir,SPMthresStruct);
%       [MapExtractStruct] = GetParamsFromSPMmat(); %select everything manually
%       [MapExtractStruct] = GetParamsFromSPMmat(SPMmatDir); %select contrast and statistics settings manually
%
%V2.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V2.0: (27.04.2015): COMPLETE OVERHAUL OF V1.0!!! (addons taken from C.Roppelt's tool "save_ThresholdedSPM") V1.0: (30.12.2014): initial implementation

%% inputs?
UseResorting = 0; %0==do not resort coordinates. (NB: doesn't change anything about clustering, but distance matrix is easier to "understand".)

%% check path to SPM.mat
try
    if(iscellstr(SPMmatDir))
        SPMmatDir = checkSPMmatDir(SPMmatDir{1});
    else
        if(ischar(SPMmatDir)) %might be a path?
            SPMmatDir = checkSPMmatDir(SPMmatDir);
        else
            %? --> select
            SPMmatDir = checkSPMmatDir([]);
        end
    end
catch CATCH_SPMmatDir
    %? --> select
    disp_catch(CATCH_SPMmatDir,[mfilename,'>checkSPMmatDir'],'CATCH_SPMmatDir')
    SPMmatDir = checkSPMmatDir([]);
end

%% get statistics via xSPM for a contrast from SPM.mat 
%% set xSPM struct 
try
    xSPM = set_xSPMstruct(SPMmatDir,SPMthresStruct.Contrast,SPMthresStruct.MultTestCor,SPMthresStruct.p,SPMthresStruct.k);
catch CATCH_set_xSPMstruct
    %everything is input by user
    disp_catch(CATCH_set_xSPMstruct,[mfilename,'>set_xSPMstruct'],'CATCH_set_xSPMstruct')
    xSPM = struct('Im', []);
    xSPM.swd   = SPMmatDir;
end

%% "update" xSPM struct
% clear SPM's "Interavtive" window as this might be needed by 'spm_getSPM'
spm_clf('Interactive');
CurrDir = pwd; %store current dir
[SPM,xSPM] = spm_getSPM(xSPM); %"update" xSPM struct
cd(CurrDir); %return to current dir

%% assign values
MapExtractStruct.Positive.Coords_mm = xSPM.XYZmm';
MapExtractStruct.Positive.StatsVals = xSPM.Z';
MapExtractStruct.Positive.Coords_vox= xSPM.XYZ';
MapExtractStruct.Positive.ConNr     = xSPM.Ic;

MapExtractStruct.Negative = [];
% MapExtractStruct.Negative.Coords_mm = [];
% MapExtractStruct.Negative.StatsVals = [];
% MapExtractStruct.Negative.Coords_vox= [];

%% apply resorting
if(UseResorting)
    FinalResortingIndices = ZYXresort(MapExtractStruct.Positive.Coords_mm);
    MapExtractStruct.Positive.Coords_mm = MapExtractStruct.Positive.Coords_mm( FinalResortingIndices,:);
    MapExtractStruct.Positive.StatsVals = MapExtractStruct.Positive.StatsVals( FinalResortingIndices);
    MapExtractStruct.Positive.Coords_vox= MapExtractStruct.Positive.Coords_vox(FinalResortingIndices,:);
end

%% write remaining info to ouput
%volume info
MapExtractStruct.V_map       = xSPM.Vspm;
MapExtractStruct.V_map.fname = [SPMmatDir,filesep,MapExtractStruct.V_map.fname];
%thresholds
MapExtractStruct.Thresholds  = [nan,xSPM.u]; %no negative only positive

end


%% subfunctions

%% checkSPMmatDir
function SPMmatDir = checkSPMmatDir(SPMmatDir)
% This function checks what is in the variable SPMmatDir, i.e. directory or filepath and does the correct assignment.
if(isempty(SPMmatDir))
    SPMmatDir = spm_select(1,'SPM.mat','Select SPM.mat for extraction of significant values from Statistics-Map...');
    SPMmatDir = fileparts(SPMmatDir);
else
    if(~exist(SPMmatDir,'dir'))
        if(exist(SPMmatDir,'file')==2) %might be path to SPM.mat instead of SPM dir
            [SPMmatDir_tmp,fname,ext] = fileparts(SPMmatDir);
            if(strcmp(fname,'SPM')&&strcmp(ext,'.mat')) %correct/as expected --> assign
                SPMmatDir = SPMmatDir_tmp;
            else
                %not correct/as expected --> select
                SPMmatDir = spm_select(1,'SPM.mat','Select SPM.mat for extraction of significant values from Statistics-Map...');
                SPMmatDir = fileparts(SPMmatDir);
            end
        else
            %not a file or directory??? --> select
            SPMmatDir = spm_select(1,'SPM.mat','Select SPM.mat for extraction of significant values from Statistics-Map...');
            SPMmatDir = fileparts(SPMmatDir);
        end
    else
        %it is a directory --> check that SPM.mat is there.
        if(~exist([SPMmatDir,filesep,'SPM.mat'],'file')) %check that SPM.mat is in the directory.
            SPMmatDir = spm_select(1,'SPM.mat','Select SPM.mat for extraction of significant values from Statistics-Map...');
            SPMmatDir = fileparts(SPMmatDir);
        else
            return;
        end
    end
end

end


%% set_xSPMstruct
function xSPM = set_xSPMstruct(SPMmatDir,Contrast,MultTestCor,p,k)
% This function adds all necessary information to the xSPM struct such that automatic selection of 
% contrast and statistics settings is done, otherwise the options will be asked to be input.

%% init 'xSPM'
xSPM = struct('Im', []);
xSPM.swd   = SPMmatDir;

%% load SPM.mat & check it
load([SPMmatDir,filesep,'SPM.mat']);

% check that our 'SPM' is a structure, that contains another structure
% called 'xCon', that contains a field called 'name'. If this is not the
% case, something's seriously wrong with that SPM structure!
if ~isstruct(SPM)
    error('Object "SPM" (loaded from file %s) is not a structure.', [SPMmatDir,filesep,'SPM.mat'])
elseif ~isfield(SPM, 'xCon')
    error('Object "SPM" does not contain required field "xCon".')
elseif ~isstruct(SPM.xCon)
    error('Field "SPM.xCon" must be a structure (but is not).')
elseif ~isfield(SPM.xCon, 'name')
    error('Expected field "name" not available in "SPM.xCon".')
end

%% find contrast
% find index and name referring to 'contrast', where 'contrast' may either
% be a string (interpreted as the contrast's name; or a single numeric
% value (interpreted as the contrast's index in 'SPM.xCon').
if(ischar(Contrast))
    matches = []; %init search empty
    for idx = 1:numel(SPM.xCon)
        if(strcmp(SPM.xCon(idx).name, Contrast)) %require perfect match!
            matches(end+1) = idx;
        end
    end
    
    if(numel(matches) == 1) %require unique perfect match!
        cidx  = matches(1);
        cname = SPM.xCon(cidx).name;
    elseif(isempty(matches))
        error('No Contrast named "%s" found.', Contrast)
    else
        error('Found more than one Contrast named "%s".', Contrast)
    end
elseif(isnumeric(Contrast) && numel(Contrast)==1)
    if(Contrast < 1)
        error('"Contrast" is smaller than 1 and thus not a valid index.')
    elseif(Contrast > numel(SPM.xCon))
        error('"Contrast" is larger than number of elements in "SPM.xCon".')
    else
        cidx  = Contrast;
        cname = SPM.xCon(cidx).name;
    end
else
    error('"Contrast" must either be a string or a single numeric value.')
end
xSPM.title = cname;
xSPM.Ic    = cidx;

%% which multiple correction and which p & k???
if(~strcmpi(MultTestCor,'FWE')&&~strcmpi(MultTestCor,'FDR')&&~strcmpi(MultTestCor,'NONE'))
    error(['Multiple Correction Type "',MultTestCor,'" is unknown!']);
else
    if(strcmpi(MultTestCor,'FDR')) %allow user to continue but give a warning text
        disp('WARNING: FDR might not be supported by SPM if the spm_defaults.m file has not been changed to allow voxel-wise FDR.');
    end
end
xSPM.thresDesc = MultTestCor;

if(isempty(p)&&~isempty(MultTestCor))
    switch(MultTestCor)
        case {'FWE','fwe'}
            p = 0.05;
        case {'FDR','fdr'}
            p = 0.01;
        case {'NONE','none'}
            p = 0.001;
        otherwise
            p = []; %set manually
    end
end
xSPM.u = p;

if(~isempty(k))
    if(mod(k,1)~=0)
        k = round(k);
    end
end
xSPM.k = k;

end