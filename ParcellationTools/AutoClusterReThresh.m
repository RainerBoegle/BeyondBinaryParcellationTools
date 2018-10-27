function [ReThreshMapPaths,ClusterMapPaths,SettingsStruct] = AutoClusterReThresh(ThreshMapPaths)
% This function allows to rethreshold maps with little ease.
%
%Usage:
%       [ReThreshMapPaths,ClusterMapPaths,SettingsStruct] = AutoClusterReThresh(ThreshMapPaths);
%       [ReThreshMapPaths,ClusterMapPaths,SettingsStruct] = AutoClusterReThresh(); %manual selection
%       [ReThreshMapPaths,ClusterMapPaths,SettingsStruct] = AutoClusterReThresh(pwd); %manual selection starting in current directory
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (20.02.2016): initial implementation based on previous tests

%% ask user to do selection
if(~exist('ThreshMapPaths','var'))
    ThreshMapPaths = cellstr(spm_select(Inf,'image','Select Maps for rethresholding...'));
else
    if(~iscellstr(ThreshMapPaths))
        if(~ischar(ThreshMapPaths))
            error('"ThreshMapPaths" must either be a cellstring or char/string!');
        else
            if(isdir(ThreshMapPaths))
                ThreshMapPaths = cellstr(spm_select(Inf,'image','Select Maps for rethresholding...',{},ThreshMapPaths));
            else
                ThreshMapPaths = cellstr(ThreshMapPaths);
            end
        end
    end
end

%% Search Distance
answerSearchDist_mm = inputdlg({'SearchDistance [mm](>0)'},'Search Distance',1,{'8'});
SearchDist_mm = str2num(answerSearchDist_mm{1}); %mm
if(SearchDist_mm<=0)
    error('SearchDist_mm must be greater zero! (Distance in mm)');
end

%% Use Connectivity Matrix
ChoiceConnectMat = questdlg('Use Connectivity Matrix to cluster connected voxels?','Use ConnectMat?','Yes','No','No');
switch(ChoiceConnectMat)
    case 'Yes'
        UseConnectMat = 1; %use connect matrix
    otherwise
        UseConnectMat = 0; %Do not use connect matrix
end

%% Retrain percentage?
answerRetainPercentage = inputdlg({'RetainPercentage [10%==0.1](>0&<1)'},'RetainPercentage',1,{'0.1'});
RetainPercentage = str2num(answerRetainPercentage{1}); %percent
if(RetainPercentage<=0)
    error('RetainPercentage must be greater zero (& smaller than 1)! (Percentage of top voxels to keep.)');
else
    if(RetainPercentage>=1)
        if(RetainPercentage>=5)
            disp('WARNING: Will assume that percentage has be entered --> dividing by 100.');
            RetainPercentage = RetainPercentage/100;
        else
            disp('WARNING: All voxels are retained after clustering!?');
        end
    end        
end

%% RetainMode? 
RetainMode = questdlg({'Which retain mode do you want to use?'; 'Use "voxels" percentage i.e. the top X% of voxels or X% of Peak Value as threshold?'},'Retain Mode?','Voxels','Threshold','Voxels');

%% XtraInputs = {}; || XtraInputs = {{'thresh',[4,-3]},{'CLsize',18}};
XtraInputs = {}; %init

%% Threshold Maps initially?
ChoiceInitThresh = questdlg('Apply initial thresholding to input map?','Initial thresholding?','Yes','No','No');
switch(ChoiceInitThresh)
    case 'Yes'
        UseInitThresh = 1; %use InitThresh
    otherwise
        UseInitThresh = 0; %Do not use InitThresh
end
if(UseInitThresh)
    answerInitThresh = inputdlg({'Initial Thresholds [Positive,Negative]'},'RetainPercentage',1,{'[4,-3]'});
    InitThresh = str2num(answerInitThresh{1}); %percent
    if(InitThresh(1)<=0)
        error('Initial Thresholds for positive side is given as negative! (Must be positive, -as the name says!)');
    end
    if(InitThresh(2)>=0)
        error('Initial Thresholds for negative side is given as positive! (Must be negative, -as the name says!)');
    end
    
    XtraInputs = {'thresh',[InitThresh(1),InitThresh(2)]};
end

%% Cluster size
ChoiceCLsize = questdlg('Use Cluster Size threshold for throwing out clusters?','Use Cluster Size Threshold?','Yes','No','No');
switch(ChoiceCLsize)
    case 'Yes'
        UseCLsizeThresh = 1; %use CLsizeThresh
    otherwise
        UseCLsizeThresh = 0; %Do not use CLsizeThresh
end
if(UseCLsizeThresh)
    answerCLsizeThresh = inputdlg({'ClusterSize Thresholds (>0)'},'RetainPercentage',1,{'18'});
    CLsizeThresh = str2num(answerCLsizeThresh{1}); %percent
    if(CLsizeThresh<=0)
        error('ClusterSize Thresholds must be positive!');
    end
    if(isempty(XtraInputs))
        XtraInputs    = {'CLsize',CLsizeThresh};
    else
        XtraInputs{2} = {'CLsize',CLsizeThresh};
    end
end

%% apply to inputs
ReThreshMapPaths= cell(length(ThreshMapPaths),1);
ClusterMapPaths = cell(length(ThreshMapPaths),1);
for Ind = 1:length(ThreshMapPaths)
    if(~isempty(XtraInputs))
        [ReThreshMapPaths{Ind},ClusterMapPaths{Ind}] = ClusterReThresh(ThreshMapPaths{Ind},SearchDist_mm,UseConnectMat,RetainPercentage,RetainMode,XtraInputs);
    else
        [ReThreshMapPaths{Ind},ClusterMapPaths{Ind}] = ClusterReThresh(ThreshMapPaths{Ind},SearchDist_mm,UseConnectMat,RetainPercentage,RetainMode);
    end
end

%% save settings structure
SettingsStruct.ThreshMapPaths  = ThreshMapPaths;
SettingsStruct.SearchDist_mm   = SearchDist_mm;
SettingsStruct.UseConnectMat   = UseConnectMat;
SettingsStruct.RetainPercentage= RetainPercentage;
SettingsStruct.RetainMode      = RetainMode;
SettingsStruct.XtraInputs      = XtraInputs;

%% Done.
disp('Done with ReThresholding by rethreholding clusters.');

end