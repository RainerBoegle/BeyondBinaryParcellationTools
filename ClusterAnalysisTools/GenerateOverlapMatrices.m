function OverlapStruct = GenerateOverlapMatrices(ClusterAnalysisDir,ControlStr)
% This function will generate the Overlap Matrices for a given ClusterAnalysis Directory,
% i.e., using the 4D-NIFTI files "ThreshICs.nii" and "Clustered<SearchDistance>mmThreshICs.nii".
%
% The Overlap Matrices include
% 1.Dice Coefficient 
%               "DiceCoeff" = 2.*sum(MaskCl_1.*MaskCl_2)./(sum(MaskCl_1)+sum(MaskCl_2))
% 2.Minimum Cluster Overlap 
%               "MinClOverlap" = sum(MaskCl_1.*MaskCl_2)./min(sum(MaskCl_1),sum(MaskCl_2))
% 3.Overlap Size (diagonal is cluster size in voxels) 
%               "OverlapSize" = sum(MaskCl_1.*MaskCl_2)
% 4.Jaccard Coefficient (similar to dice coefficient but considers overlapping areas differently in numberator and denominator normalization.)
%               "JaccardCoeff" = sum(MaskCl_1.*MaskCl_2)./(sum(MaskCl_1)+sum(MaskCl_2)-sum(MaskCl_1.*MaskCl_2))
% 5.Sum Normalized Statistic Overlap in Cluster Mask, i.e., for a cluster all statistics values
%   are summed and each value divided by the sum, normalizing the sum to 1.
%   Therefore this overlap with another cluster, is the amount of statistics values that are
%   overlapping the other cluster while including the (normalized) significance amplitude! 
%   In contradistinction to the DiceCoeff which does not include significance amplitude!
%   NB: The other cluster is taken as a mask, i.e. its values are all set to 1!!! 
%       Not both taken as statistics values! 
%       --> Therefore this matrix is NOT SYMMETRIC!!!
%               "SumContribution" = sum(NormStatsValsCl_1.*MaskCl_2);
% 6.Absolute value of the Correlation coefficient (Pearson correlation coefficient of statistics values!)
%               "AbsCorrCoeff" = abs(corr(StatsValsCl_1,StatsValsCl_2));
%
%
%
%Usage:
%       OverlapStruct = GenerateOverlapMatrices(ClusterAnalysisDir,ControlStr);
%       OverlapStruct = GenerateOverlapMatrices(); %select ClusterAnalysis directory using spm_select AND GENERATE OverlapStruct and save it in ClusterAnalysis directory.
%       OverlapStruct = GenerateOverlapMatrices([],'gen'); %the same as above
%       OverlapStruct = GenerateOverlapMatrices(ClusterAnalysisDir,'load'); %load OverlapStruct from OverlapStruct.mat if available in ClusterAnalysis directory.
%       OverlapStruct = GenerateOverlapMatrices([],'load'); %load OverlapStruct from OverlapStruct.mat if available in ClusterAnalysis directory.
%
%
%V1.1
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.1: (05.11.2016): add Jaccard coefficient V1.0: (02.11.2016): initial implementation based on previous tests

%% Check inputs
%ClusterAnalysisDir
if(~exist('ClusterAnalysisDir','var'))
    ClusterAnalysisDir = spm_select(1,'dir','Select ClusterAnalysis Directory...');
else
    if(isempty(ClusterAnalysisDir))
        ClusterAnalysisDir = spm_select(1,'dir','Select ClusterAnalysis Directory...');
    else
        if(~ischar(ClusterAnalysisDir))
            if(iscellstr(ClusterAnalysisDir))
                ClusterAnalysisDir = ClusterAnalysisDir{1}; %take first
            else
                error('ClusterAnalysisDir must be a string that points to the ClusterAnalysis directory, i.e. the directory that contains the NIFTI-files "ThreshICs.nii" and "Clustered<SearchDistance>mmThreshICs.nii".');
            end
        end
    end
end

%ControlStr
if(~exist('ControlStr','var'))
    ControlStr = 'gen';
else
    if(isempty(ControlStr))
        ControlStr = 'gen';
    else
        if(~ischar(ControlStr))
            error('ControlStr must be a string! (Either "gen" in case the Overlaps should be generated and "load" in case it should be loaded from the "OverlapStruct.mat" in the ClusterAnalysis Directory (if available).');
        end
    end
end
switch(ControlStr)
    case {'gen','Gen','Generate','generate','g','G'}
        disp(['Will generate overlap matrices and save "OverlapStruct.mat" in ClusterAnalysis directory "',ClusterAnalysisDir,'".']);
    case {'load','Load','l','L'}
        disp(['Will load overlap matrices by loading "OverlapStruct.mat" from ClusterAnalysis directory "',ClusterAnalysisDir,'".']);
        if(exist([ClusterAnalysisDir,filesep,'OverlapStruct.mat'],'file'))
            load([ClusterAnalysisDir,filesep,'OverlapStruct.mat']);
            disp(' ');
            disp('Done.');
            disp(' ');
            return;
        else
            error(['Could NOT FIND "OverlapStruct.mat" in directory "',ClusterAnalysisDir,'"']);
        end
    otherwise
        error(['Unknown ControlStr "',ControlStr,'"!']);
end

%% get NIFTI-files from ClusterAnalysisDir
ThreshICsPath          = cellstr(spm_select('FPList',ClusterAnalysisDir,'^ThreshICs.nii'));
ClusteredThreshICsPath = cellstr(spm_select('FPList',ClusterAnalysisDir,'^Clustered.*ThreshICs.nii'));

%ThreshICsPath
if(length(ThreshICsPath)>1)
    [SelInd,ok] = listdlg('ListString',ThreshICsPath,'SelectionMode','single','Name','Select Stats-file','PromptString','Select Statistic file for analysis...');
    if(~ok)
        disp('Quit by user choice.');
        OutputDir = [];
        ClusterEncodeStruct = [];
        return;
    else
        ThreshICsPath = ThreshICsPath{SelInd};
    end
elseif(length(ThreshICsPath)==1)
    ThreshICsPath = ThreshICsPath{1};
else
    error(['Could not find statistics file in directory "',ClusterAnalysisDir,'"...']);
end

%ClusteredThreshICsPath
if(length(ClusteredThreshICsPath)>1)
    [SelInd,ok] = listdlg('ListString',ClusteredThreshICsPath,'SelectionMode','single','Name','Select Cluster-file','PromptString','Select Clustering file for analysis...');
    if(~ok)
        disp('Quit by user choice.');
        OutputDir = [];
        ClusterEncodeStruct = [];
        return;
    else
        ClusteredThreshICsPath = ClusteredThreshICsPath{SelInd};
    end
elseif(length(ClusteredThreshICsPath)==1)
    ClusteredThreshICsPath = ClusteredThreshICsPath{1};
else
    error(['Could not find Clustering file in directory "',ClusterAnalysisDir,'"...']);
end

%% get 4D Data from files
disp(['Loading data from directory "',ClusterAnalysisDir,'"...']);
V_ThreshICs          = spm_vol(ThreshICsPath);
V_ClusteredThreshICs = spm_vol(ClusteredThreshICsPath);
if(length(V_ClusteredThreshICs)==length(V_ThreshICs)&&length(V_ThreshICs)~=1)
    NVolumes = length(V_ThreshICs);
end

%ThreshICs4D          = V_ThreshICs(1).private.dat(:,:,:,:); %the first one is fine enough, because we can access all data from there.
%ClusteredThreshICs4D = V_ClusteredThreshICs(1).private.dat(:,:,:,:); %the first one is fine enough, because we can access all data from there.
ThreshICs4D          = zeros(V_ThreshICs(1).dim(1),V_ThreshICs(1).dim(2),V_ThreshICs(1).dim(3),NVolumes);
ClusteredThreshICs4D = zeros(V_ClusteredThreshICs(1).dim(1),V_ClusteredThreshICs(1).dim(2),V_ClusteredThreshICs(1).dim(3),NVolumes);
for IndVol = 1:NVolumes
    ThreshICs4D(:,:,:,IndVol)          = V_ThreshICs(IndVol).private.dat(:,:,:,V_ThreshICs(IndVol).n(1));
    ClusteredThreshICs4D(:,:,:,IndVol) = V_ClusteredThreshICs(IndVol).private.dat(:,:,:,V_ClusteredThreshICs(IndVol).n(1));
end
if(any(ClusteredThreshICs4D(:)<0))
    error('ClusteredThreshICs4D, i.e., all the clusters in space per IC contain negative numbers! This can not be, there must have been something wrong with the generation of the "Clustered<SearchDistance>mmThreshICs.nii" file.');
end

%% create mask and get NVoxels (in mask)
disp('Generating Mask...');
Mask3D = (sum(ClusteredThreshICs4D,4)>0 | sum(abs(ThreshICs4D),4)>0);
NVoxels= length(find(Mask3D(:)~=0));

%% get number of IC "NICs" and use mask to make 4D Data into 2D Voxels-x-ICs data matrices 
disp('Extract data using Mask...');
NICs = size(ClusteredThreshICs4D,4);
ThreshICs2D         = zeros(NVoxels,NICs);
ClusteredThreshICs2D= zeros(NVoxels,NICs);
for IndIC = 1:NICs
    CurrThreshICs      = ThreshICs4D(:,:,:,IndIC);
    ClusteredThreshICs = ClusteredThreshICs4D(:,:,:,IndIC);
    ThreshICs2D(:,IndIC)         = CurrThreshICs(     Mask3D(:)~=0);
    ClusteredThreshICs2D(:,IndIC)= ClusteredThreshICs(Mask3D(:)~=0);
    if(all(ThreshICs2D(:,IndIC)==0)||all(ClusteredThreshICs2D(:,IndIC)==0))
        error(['IC ',num2str(IndIC),' only contains zeros!']);
    end
end

%% get number of clusters per IC "ClsPerIC" and total number of clusters "totalClNum"
%% AND generate the listing of [ICnum,ClusterNum] for 1:totalClNum, i.e. the matrix "ICnumClNum" of size totalClNum-x-2 
disp('Generate listing of [ICnum,ClusterNum] for all clusters...');
ClsPerIC   = max(ClusteredThreshICs2D,[],1)';
totalClNum = sum(ClsPerIC);
ICnumClNum = zeros(totalClNum,2);
Index      = 0; %init
for IndIC = 1:NICs
    for IndCl = 1:ClsPerIC(IndIC)
        Index = Index + 1; %add this one
        ICnumClNum(Index,:) = [IndIC, IndCl];
    end
end
if(Index~=totalClNum)
    error('Something is wrong with the total number of clusters and the number of clusters per IC, -they do not add up.');
end

%% Expand 2D data matrices into 2D Voxels-x-totalClNum NB: clusters are made into masks
disp('Expand data into individual clusters...');
StatsVals2D    = zeros(NVoxels,totalClNum);
ClusterMasks2D = zeros(NVoxels,totalClNum);
for IndCl = 1:totalClNum
    ClusterMasks2D(:,IndCl) = double(ClusteredThreshICs2D(:,ICnumClNum(IndCl,1))==ICnumClNum(IndCl,2));
    StatsVals2D(:,IndCl)    = ThreshICs2D(:,ICnumClNum(IndCl,1)).*ClusterMasks2D(:,IndCl); %take only stats vals from the cluster
end

%% Normalize StatsVals
disp('Normalize stats vals via total sum == 1...');
AbsStatsVals2D    = abs(StatsVals2D);
SumAbsStatsVals2D = sum(AbsStatsVals2D,1);
NormStatsVals2D   = AbsStatsVals2D./repmat(SumAbsStatsVals2D,NVoxels,1);

%% Generate Overlap Matrices via Matrix products
disp('Generating Overlap Matrices...');
ClusterSizes1   = repmat(sum(ClusterMasks2D,1),size(ClusterMasks2D,2),1);
ClusterSizes2   = repmat(sum(ClusterMasks2D,1)',1,size(ClusterMasks2D,2));

OverlapSize     =    ClusterMasks2D' * ClusterMasks2D;
DiceCoeff       = 2.*OverlapSize./(   ClusterSizes1+ClusterSizes2);
MinClOverlap    =    OverlapSize./min(ClusterSizes1,ClusterSizes2);
JaccardCoeff    =    OverlapSize./(   ClusterSizes1+ClusterSizes2-OverlapSize);
SumContribution =    NormStatsVals2D'* ClusterMasks2D;
AbsCorrCoeff    = abs(corr(StatsVals2D,StatsVals2D));

%% assign matrices to OverlapStruct
%the settings
Setup    = struct('ClusterAnalysisDir',ClusterAnalysisDir,'ThreshICsPath',ThreshICsPath,'ClusteredThreshICsPath',ClusteredThreshICsPath,...
                  'NVoxels',NVoxels,'NICs',NICs,'totalClNum',totalClNum,'ClsPerIC',ClsPerIC,'ICnumClNum',ICnumClNum);
%the Input data
InputData= struct('Mask3D',Mask3D,'StatsVals2D',StatsVals2D,'NormStatsVals2D',NormStatsVals2D,'ClusterMasks2D',ClusterMasks2D);

%finally the OverlapStruct with the Overlap Matrices, the Setup & InputData
OverlapStruct = struct('DiceCoeff',DiceCoeff,'MinClOverlap',MinClOverlap,'OverlapSize',OverlapSize,'JaccardCoeff',JaccardCoeff,'SumContribution',SumContribution,'AbsCorrCoeff',AbsCorrCoeff,...
                       'Setup',Setup,'InputData',InputData);
                   
%% save OverlapStruct in ClusterAnalysisDir
disp(['Saving "OverlapStruct.mat" in directory "',ClusterAnalysisDir,'".']);
save([ClusterAnalysisDir,filesep,'OverlapStruct.mat'],'OverlapStruct');

%% done.
disp(' ');
disp('Done.');
disp(' ');

end
