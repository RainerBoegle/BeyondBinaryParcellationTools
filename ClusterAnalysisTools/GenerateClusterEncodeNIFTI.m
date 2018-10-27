function [OutputDir,ClusterEncodeStruct] = GenerateClusterEncodeNIFTI(ClusterAnalysisDir,OverlapMatThresholds,AnalysisNAME)
% This function can generate an encoding of the clusters as indicated by the
% Overlap Matrices stored in the "OverlapStruct" (OverlapStruct.mat) in the 
% ClusterAnalysis directory "ClusterAnalysisDir", AND the thresholds indicated
% in the structure "OverlapMatThresholds".
% See GenerateOverlapMatrices.m for details on the overlap measures, their definition and meaning.
%
% The thresholds structure can have the following fields, if any one is missing,
% it will be set to a default value indicated below.
% For any overlap to be accepted, it has to be equal or higher than this threshold.
%
%OverlapMatThresholds.
%                    .NVoxelMin            ( 10)   Minimum number of voxels per cluster & cluster overlap 
%                    .MinDiceCoeff         (0.2)   Minimum Dice coefficient per cluster overlap 
%                    .MinMinClOverlap      (0.5)   Minimum "minimum cluster overlap" per cluster overlap 
%                    .MinJaccardCoeff      (0.1)   Minimum Jaccard Coefficient per cluster overlap 
%                    .MinSumContribution   (0.5)   Minimum Sum Normalized Statistics Contribution per cluster overlap 
%                    .MinAbsCorrCoeff      (0.3)   Minimum Correlation coefficient (Pearson) of statistics values.
%
% 
%
%NB:  The "OverlapStruct" and the "OverlapMatThresholds" will be stored in
%     the "ClusterEncodeStruct" which will be saved in the Output directory
%     together with the generated NIFTI-files.
%NB2: The OutputDir will be generated from the ClusterAnalysis directory
%     and the AnalysisNAME, if AnalysisNAME is not given it will be replaced by
%     a datestring.
%
%
%
%Usage:
%       [OutputDir,ClusterEncodeStruct] = GenerateClusterEncodeNIFTI(ClusterAnalysisDir,OverlapMatThresholds,AnalysisNAME);
%
%
%
%
%V1.1
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.1: (06.11.2016): added Jaccard Coeff and AbsCorrCoeff. V1.0: (05.11.2016): initial implementation based on previous tests

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

%OverlapMatThresholds
if(~exist('OverlapMatThresholds','var'))
    disp('Will use defaults for OverlapMatThresholds...');
    OverlapMatThresholds.NVoxelMin          =  10; %.NVoxelMin            ( 10)   Minimum number of voxels per cluster & cluster overlap 
    OverlapMatThresholds.MinDiceCoeff       = 0.2; %.MinDiceCoeff         (0.2)   Minimum Dice coefficient per cluster overlap 
    OverlapMatThresholds.MinMinClOverlap    = 0.5; %.MinMinClOverlap      (0.5)   Minimum "minimum cluster overlap" per cluster overlap 
    OverlapMatThresholds.MinJaccardCoeff    = 0.1; %.MinJaccardCoeff      (0.1)   Minimum Jaccard Coefficient per cluster overlap 
    OverlapMatThresholds.MinSumContribution = 0.5; %.MinSumContribution   (0.5)   Minimum Sum Normalized Statistics Contribution per cluster overlap 
    OverlapMatThresholds.MinAbsCorrCoeff    = 0.3; %.MinAbsCorrCoeff      (0.3)   Minimum Correlation coefficient (Pearson) of statistics values.
    
    disp('Will use defaults for OverlapMatThresholds.NVoxelMin = 10;');
    disp('Will use defaults for OverlapMatThresholds.MinDiceCoeff = 0.2;');
    disp('Will use defaults for OverlapMatThresholds.MinMinClOverlap = 0.5;');
    disp('Will use defaults for OverlapMatThresholds.MinJaccardCoeff = 0.1;');
    disp('Will use defaults for OverlapMatThresholds.MinSumContribution = 0.5;');
    disp('Will use defaults for OverlapMatThresholds.MinAbsCorrCoeff = 0.3;');
else
    if(isempty(OverlapMatThresholds))
        disp('Will use defaults for OverlapMatThresholds...');
        OverlapMatThresholds.NVoxelMin          =  10; %.NVoxelMin            ( 10)   Minimum number of voxels per cluster & cluster overlap
        OverlapMatThresholds.MinDiceCoeff       = 0.2; %.MinDiceCoeff         (0.2)   Minimum Dice coefficient per cluster overlap
        OverlapMatThresholds.MinMinClOverlap    = 0.5; %.MinMinClOverlap      (0.5)   Minimum "minimum cluster overlap" per cluster overlap
        OverlapMatThresholds.MinJaccardCoeff    = 0.1; %.MinJaccardCoeff      (0.1)   Minimum Jaccard Coefficient per cluster overlap
        OverlapMatThresholds.MinSumContribution = 0.5; %.MinSumContribution   (0.5)   Minimum Sum Normalized Statistics Contribution per cluster overlap
        OverlapMatThresholds.MinAbsCorrCoeff    = 0.3; %.MinAbsCorrCoeff      (0.3)   Minimum Correlation coefficient (Pearson) of statistics values.
        
        disp('Will use defaults for OverlapMatThresholds.NVoxelMin = 10;');
        disp('Will use defaults for OverlapMatThresholds.MinDiceCoeff = 0.2;');
        disp('Will use defaults for OverlapMatThresholds.MinMinClOverlap = 0.5;');
        disp('Will use defaults for OverlapMatThresholds.MinJaccardCoeff = 0.1;');
        disp('Will use defaults for OverlapMatThresholds.MinSumContribution = 0.5;');
        disp('Will use defaults for OverlapMatThresholds.MinAbsCorrCoeff = 0.3;');
    else
        if(isstruct(OverlapMatThresholds))
            %NVoxelMin 
            if(~isfield(OverlapMatThresholds,'NVoxelMin'))
                disp('Will use defaults for OverlapMatThresholds.NVoxelMin = 10;');
                OverlapMatThresholds.NVoxelMin          =  10; %.NVoxelMin            ( 10)   Minimum number of voxels per cluster & cluster overlap
            elseif(isempty(OverlapMatThresholds.NVoxelMin))
                disp('Will use defaults for OverlapMatThresholds.NVoxelMin = 10;');
                OverlapMatThresholds.NVoxelMin          =  10; %.NVoxelMin            ( 10)   Minimum number of voxels per cluster & cluster overlap
            elseif(OverlapMatThresholds.NVoxelMin>=0)
                disp(['Will use OverlapMatThresholds.NVoxelMin == ',num2str(OverlapMatThresholds.NVoxelMin)]);
            else
                error('OverlapMatThresholds.NVoxelMin must be greater or equal to 0.');
            end
            %MinDiceCoeff
            if(~isfield(OverlapMatThresholds,'MinDiceCoeff'))
                disp('Will use defaults for OverlapMatThresholds.MinDiceCoeff = 0.2;');
                OverlapMatThresholds.MinDiceCoeff       = 0.2; %.MinDiceCoeff         (0.2)   Minimum Dice coefficient per cluster overlap
            elseif(isempty(OverlapMatThresholds.MinDiceCoeff))
                disp('Will use defaults for OverlapMatThresholds.MinDiceCoeff = 0.2;');
                OverlapMatThresholds.MinDiceCoeff       = 0.2; %.MinDiceCoeff         (0.2)   Minimum Dice coefficient per cluster overlap
            elseif(OverlapMatThresholds.MinDiceCoeff<=1&&OverlapMatThresholds.MinDiceCoeff>=0)
                disp(['Will use OverlapMatThresholds.MinDiceCoeff == ',num2str(OverlapMatThresholds.MinDiceCoeff)]);
            else
                error('OverlapMatThresholds.MinDiceCoeff must be between 0 and 1.');
            end
            %MinMinClOverlap
            if(~isfield(OverlapMatThresholds,'MinMinClOverlap'))
                disp('Will use defaults for OverlapMatThresholds.MinMinClOverlap = 0.5;');
                OverlapMatThresholds.MinMinClOverlap    = 0.5; %.MinMinClOverlap      (0.5)   Minimum "minimum cluster overlap" per cluster overlap
            elseif(isempty(OverlapMatThresholds.MinMinClOverlap))
                disp('Will use defaults for OverlapMatThresholds.MinMinClOverlap = 0.5;');
                OverlapMatThresholds.MinMinClOverlap    = 0.5; %.MinMinClOverlap      (0.5)   Minimum "minimum cluster overlap" per cluster overlap
            elseif(OverlapMatThresholds.MinMinClOverlap<=1&&OverlapMatThresholds.MinMinClOverlap>=0)
                disp(['Will use OverlapMatThresholds.MinMinClOverlap == ',num2str(OverlapMatThresholds.MinMinClOverlap)]);
            else
                error('OverlapMatThresholds.MinMinClOverlap must be between 0 and 1.');
            end
            %MinJaccardCoeff
            if(~isfield(OverlapMatThresholds,'MinJaccardCoeff'))
                disp('Will use defaults for OverlapMatThresholds.MinJaccardCoeff = 0.1;');
                OverlapMatThresholds.MinJaccardCoeff    = 0.1; %.MinJaccardCoeff      (0.1)   Minimum Jaccard Coefficient per cluster overlap
            elseif(isempty(OverlapMatThresholds.MinJaccardCoeff))
                disp('Will use defaults for OverlapMatThresholds.MinJaccardCoeff = 0.1;');
                OverlapMatThresholds.MinJaccardCoeff    = 0.1; %.MinJaccardCoeff      (0.1)   Minimum Jaccard Coefficient per cluster overlap
            elseif(OverlapMatThresholds.MinJaccardCoeff<=1&&OverlapMatThresholds.MinJaccardCoeff>=0)
                disp(['Will use OverlapMatThresholds.MinJaccardCoeff == ',num2str(OverlapMatThresholds.MinJaccardCoeff)]);
            else
                error('OverlapMatThresholds.MinJaccardCoeff must be between 0 and 1.');
            end
            %MinSumContribution
            if(~isfield(OverlapMatThresholds,'MinSumContribution'))
                disp('Will use defaults for OverlapMatThresholds.MinSumContribution = 0.5;');
                OverlapMatThresholds.MinSumContribution = 0.5; %.MinSumContribution   (0.5)   Minimum Sum Normalized Statistics Contribution per cluster overlap
            elseif(isempty(OverlapMatThresholds.MinSumContribution))
                disp('Will use defaults for OverlapMatThresholds.MinSumContribution = 0.5;');
                OverlapMatThresholds.MinSumContribution = 0.5; %.MinSumContribution   (0.5)   Minimum Sum Normalized Statistics Contribution per cluster overlap
            elseif(OverlapMatThresholds.MinSumContribution<=1&&OverlapMatThresholds.MinSumContribution>=0)
                disp(['Will use OverlapMatThresholds.MinSumContribution == ',num2str(OverlapMatThresholds.MinSumContribution)]);
            else
                error('OverlapMatThresholds.MinSumContribution must be between 0 and 1.');
            end
            %MinAbsCorrCoeff
            if(~isfield(OverlapMatThresholds,'MinAbsCorrCoeff'))
                disp('Will use defaults for OverlapMatThresholds.MinAbsCorrCoeff = 0.3;');
                OverlapMatThresholds.MinAbsCorrCoeff       = 0.3; %.MinAbsCorrCoeff         (0.3)   Minimum Correlation coefficient (Pearson) of statistics values.
            elseif(isempty(OverlapMatThresholds.MinAbsCorrCoeff))
                disp('Will use defaults for OverlapMatThresholds.MinAbsCorrCoeff = 0.3;');
                OverlapMatThresholds.MinAbsCorrCoeff       = 0.3; %.MinAbsCorrCoeff         (0.3)   Minimum Correlation coefficient (Pearson) of statistics values.
            elseif(OverlapMatThresholds.MinAbsCorrCoeff<=1&&OverlapMatThresholds.MinAbsCorrCoeff>=0)
                disp(['Will use OverlapMatThresholds.MinAbsCorrCoeff == ',num2str(OverlapMatThresholds.MinAbsCorrCoeff)]);
            else
                error('OverlapMatThresholds.MinAbsCorrCoeff must be between 0 and 1.');
            end
        else
            error('OverlapMatThresholds must be struct! (Fields may include .NVoxelMin .MinDiceCoeff .MinMinClOverlap .MinJaccardCoeff .MinSumContribution .MinAbsCorrCoeff)');
        end
    end
end

%AnalysisNAME
if(~exist('AnalysisNAME','var'))
    AnalysisNAME = ['ClusterEncode',datestr(now,'yyyymmmdd_HHMM')];
else
    if(isempty(AnalysisNAME))
        AnalysisNAME = ['ClusterEncode',datestr(now,'yyyymmmdd_HHMM')];
    else
        if(~ischar(AnalysisNAME))
            if(iscellstr(AnalysisNAME))
                AnalysisNAME = AnalysisNAME{1}; %take first
            else
                error('AnalysisNAME has to be string!');
            end
        end
    end
end

%% load OverlapStruct
disp(['Will load "OverlapStruct" from "OverlapStruct.mat" in ClusterAnalysis directory "',ClusterAnalysisDir,'".']);
if(exist([ClusterAnalysisDir,filesep,'OverlapStruct.mat'],'file'))
    load([ClusterAnalysisDir,filesep,'OverlapStruct.mat']);
else
    error(['Could NOT FIND "OverlapStruct.mat" in directory "',ClusterAnalysisDir,'"']);
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
disp(['Loading NIFTI data from directory "',ClusterAnalysisDir,'"...']);
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

%% Check that this fits the OverlapStruct
if((totalClNum~=OverlapStruct.Setup.totalClNum)||(totalClNum~=size(OverlapStruct.DiceCoeff,1)))
    error('Total number of clusters as loaded from NIFTI-files does not fit the number of clusters indicated in the OverlapStruct!');
end
for IndIC = 1:NICs
    if(ClsPerIC(IndIC)~=OverlapStruct.Setup.ClsPerIC(IndIC))
        error(['IC ',num2str(IndIC),': number of clusters (',num2str(ClsPerIC(IndIC)),') does not fit the number of clusters (',num2str(OverlapStruct.Setup.ClsPerIC(IndIC)),') indicated in the OverlapStruct!']);
    end
end
for IndCL = 1:totalClNum
    if(ICnumClNum(Index,1)~=OverlapStruct.Setup.ICnumClNum(Index,1))
        error(['Cluster ',num2str(IndCL),': IC ',num2str(ICnumClNum(Index,1)),' mismatch to IC ',num2str(OverlapStruct.Setup.ICnumClNum(Index,1)),' indicated in the OverlapStruct!']);
    end
    if(ICnumClNum(Index,2)~=OverlapStruct.Setup.ICnumClNum(Index,2))
        error(['Cluster ',num2str(IndCL),': IC ',num2str(ICnumClNum(Index,1)),' SubCluster ',num2str(ICnumClNum(Index,2)),' mismatch to SubCluster ',num2str(OverlapStruct.Setup.ICnumClNum(Index,2)),' of IC ',num2str(OverlapStruct.Setup.ICnumClNum(Index,1)),' indicated in the OverlapStruct!']);
    end
end
if(size(Mask3D,1)~=size(OverlapStruct.InputData.Mask3D,1))
    error('Mask3D does not fit Mask3D indicated in the OverlapStruct!');
end
if(size(Mask3D,2)~=size(OverlapStruct.InputData.Mask3D,2))
    error('Mask3D does not fit Mask3D indicated in the OverlapStruct!');
end
if(size(Mask3D,3)~=size(OverlapStruct.InputData.Mask3D,3))
    error('Mask3D does not fit Mask3D indicated in the OverlapStruct!');
end        

%% Select cluster assignments per IC, i.e. per IC only one cluster can correspond to the cluster that is investigated.
%% The Clusters have to survive the NVoxelMin threshold AND Minimum MinClusterOverlap MinMinClOverlap threshold, 
disp('Collect overlaps of clusters...');
damn this is not enough#1#1!!!
ClusterICoverlap = zeros(totalClNum,NICs);
NVoxelClOverlap  = zeros(totalClNum,NICs);
for IndCl1 = 1:totalClNum
    [ClusterICoverlap(IndCl1,:), NVoxelClOverlap(IndCl1,:)] = findBestClusterPerIC(OverlapStruct,ICnumClNum,IndCl1,OverlapMatThresholds);
end
keyboard;    

%% apply thresholds to OverlapStruct Overlap Matrices and generate one thresholded OverlapMat
OverlapMat = (OverlapStruct.DiceCoeff>=OverlapMatThresholds.MinDiceCoeff).*...
             (OverlapStruct.JaccardCoeff>=OverlapMatThresholds.MinJaccardCoeff).*...
             (OverlapStruct.MinClOverlap>=OverlapMatThresholds.MinMinClOverlap).*...
             (OverlapStruct.SumContribution>=OverlapMatThresholds.MinSumContribution).*...
             (OverlapStruct.AbsCorrCoeff>=OverlapMatThresholds.MinAbsCorrCoeff).*...
             (OverlapStruct.OverlapSize>=OverlapMatThresholds.NVoxelMin);
         
%% Symmetrize OverlapMat?
%OverlapMat = (OverlapMat + OverlapMat')~=0;

%% Generate ClusterEncode 2D Voxel-x-ClusterEncodeNum matrix for each row in OverlapMat
disp('Treating clusters according to OverlapMat...');
AllClusterSizes = diag(OverlapStruct.OverlapSize);
EncodeIndices = cell(totalClNum,1);
ClusterSizes  = cell(totalClNum,1); %sizes of clusters that were combined correspond with EncodeIndices
NClusters     = zeros(totalClNum,1); %number of clusters per ClusterEncodeNum, this can be zero if none survived, 1 for unique or up to N
%also need significance level in case clusters do overlap minimally...damn

ClusterEncode2D = zeros(size(ClusterMasks2D));
for IndCL = 1:totalClNum
    CurrData = OverlapMat(IndCL,:);
    if(all(CurrData(:)==0)) %if all zeros continue to next row 
        disp(['Cluster ',num2str(IndCL,['%0',num2str(ceil(log10(totalClNum))),'d']),' of ',num2str(totalClNum),' is empty!']);
        continue;
    else
        if(any(CurrData(:)~=0)) %if any nonzero then find indices and encode the overlap into the volume of the row index and save the indices 
            EncodeIndices{IndCL} = find(CurrData~=0);
            [ClusterSizes{IndCL},SortInds] = sort(AllClusterSizes(EncodeIndices{IndCL}),'descend'); %sort cluster sizes descending
            EncodeIndices{IndCL} = EncodeIndices{IndCL}(SortInds); %resort by cluster size
            NClusters(IndCL)     = length(EncodeIndices{IndCL});
            disp(['Cluster ',num2str(IndCL,['%0',num2str(ceil(log10(totalClNum))),'d']),' of ',num2str(totalClNum),' contains ',num2str(NClusters(IndCL)),' Clusters of sizes [',num2str(ClusterSizes{IndCL}(:)'),'].']);
            
            ClusterEncode2D(:,IndCL) = sum(ClusterMasks2D(:,EncodeIndices{IndCL}),2);
        end
    end
end

%% Collect raw cluster results
RawClusterEncode = struct('ClusterEncode2D',ClusterEncode2D,'EncodeIndices',{EncodeIndices},'ClusterSizes',{ClusterSizes},'NClusters',NClusters);

%% Remove zeros and re-evaluate overlaps using MinClusterOverlap
%sort max overlap to min (which is zero at this point)
[NClusters,SortInds] = sort(NClusters,'descend');
EncodeIndices   = EncodeIndices(SortInds);
ClusterSizes    = ClusterSizes(SortInds);
ClusterEncode2D = ClusterEncode2D(:,SortInds);

%remove zeros i.e. thin out
EncodeIndices(    NClusters==0) = [];
ClusterSizes(     NClusters==0) = [];
ClusterEncode2D(:,NClusters==0) = [];
NClusters(        NClusters==0) = [];

% %for all overlapping ones remove those voxels that are 1 i.e. unique
% ClusterEncode2D(:,NClusters>1) = ClusterEncode2D(:,NClusters>1).*(ClusterEncode2D(:,NClusters>1)>1); %mask out uniques in overlapping ones
% %ClusterEncode2D(:,NClusters>1) = ClusterEncode2D(:,NClusters>1).*repmat(sum(ClusterEncode2D(:,NClusters==1),2)>1,1,length(find(NClusters>1))); %mask out the uniques from here.
% 
% %for all uniques remove overlapping voxels
% %ClusterEncode2D(:,NClusters==1) = ClusterEncode2D(:,NClusters==1).*repmat(sum(ClusterEncode2D(:,NClusters>1),2)<2,1,length(find(NClusters==1))); %remove those that are overlapping ones
% ClusterEncode2D(:,NClusters==1) = ClusterEncode2D(:,NClusters==1).*repmat(sum(ClusterEncode2D(:,NClusters>1),2)>1,1,length(find(NClusters==1))); %remove those that are uniques ones
% ClusterEncode2D(:,NClusters==1) = ClusterEncode2D(:,NClusters==1).*(ClusterEncode2D(:,NClusters==1)>1); %mask out uniques that are remaining

%MinClusterOverlap
ClusterSizes1   = repmat(sum(ClusterEncode2D~=0,1),size(ClusterEncode2D~=0,2),1);
ClusterSizes2   = repmat(sum(ClusterEncode2D~=0,1)',1,size(ClusterEncode2D~=0,2));

OverlapSize     = double(ClusterEncode2D~=0)' * double(ClusterEncode2D~=0);
MinClOverlap    = OverlapSize./min(ClusterSizes1,ClusterSizes2);

figure(421); clf; 
subplot(1,2,1); imagesc(MinClOverlap,[0 1]); title('MinClOverlap of thinned out Cluster Encode');
set(gca,'YTick',1:length(NClusters)); set(gca,'YTickLabel',cellstr(num2str(NClusters(:)))); set(gca,'XTick',1:length(NClusters)); set(gca,'XTickLabel',cellstr(num2str(NClusters(:)))); set(gca,'XTickLabelRotation',-45);
subplot(1,2,2); imagesc(sum(MinClOverlap>=OverlapMatThresholds.MinMinClOverlap,2)); title(['sum(MinClOverlap>=',num2str(OverlapMatThresholds.MinMinClOverlap),',2) of thinned out Cluster Encode']);
set(gca,'YTick',1:length(NClusters)); set(gca,'YTickLabel',cellstr(num2str(NClusters(:)))); 

figure(422); clf; 
subplot(1,2,1); imagesc(OverlapSize); title('OverlapSize of thinned out Cluster Encode');
set(gca,'YTick',1:length(NClusters)); set(gca,'YTickLabel',cellstr(num2str(NClusters(:)))); set(gca,'XTick',1:length(NClusters)); set(gca,'XTickLabel',cellstr(num2str(NClusters(:)))); set(gca,'XTickLabelRotation',-45);
subplot(1,2,2); imagesc(sum(OverlapSize>=OverlapMatThresholds.NVoxelMin,2)); title(['sum(OverlapSize>=',num2str(OverlapMatThresholds.NVoxelMin),',2) of thinned out Cluster Encode']);
set(gca,'YTick',1:length(NClusters)); set(gca,'YTickLabel',cellstr(num2str(NClusters(:)))); 

figure(42); clf; 
subplot(1,2,1); imagesc((MinClOverlap>=OverlapMatThresholds.MinMinClOverlap).*(OverlapSize>=OverlapMatThresholds.NVoxelMin),[0 1]); title('MinClOverlap.*OverlapSize THRESHOLDED of thinned out Cluster Encode');
set(gca,'YTick',1:length(NClusters)); set(gca,'YTickLabel',cellstr(num2str(NClusters(:)))); set(gca,'XTick',1:length(NClusters)); set(gca,'XTickLabel',cellstr(num2str(NClusters(:)))); set(gca,'XTickLabelRotation',-45);
subplot(1,2,2); imagesc(sum((MinClOverlap>=OverlapMatThresholds.MinMinClOverlap).*(OverlapSize>=OverlapMatThresholds.NVoxelMin),2)); title(['sum(OverlapSize>=',num2str(OverlapMatThresholds.NVoxelMin),' .* MinClOverlap>=',num2str(OverlapMatThresholds.MinMinClOverlap),',2) of thinned out Cluster Encode']);
set(gca,'YTick',1:length(NClusters)); set(gca,'YTickLabel',cellstr(num2str(NClusters(:)))); 


%% Collect thinned out Cluster Encode
ClusterEncode = struct('ClusterEncode2D',ClusterEncode2D,'EncodeIndices',{EncodeIndices},'ClusterSizes',{ClusterSizes},'NClusters',NClusters,'OverlapSize',OverlapSize,'MinClOverlap',MinClOverlap);

%% setup ClusterEncodeStruct
ClusterEncodeStruct = struct('OverlapMat',OverlapMat,'OverlapStruct',OverlapStruct,'OverlapMatThresholds',OverlapMatThresholds,'RawClusterEncode',RawClusterEncode,'ClusterEncode',ClusterEncode);

%% create OutputDir
OutputDir = [ClusterAnalysisDir,filesep,AnalysisNAME];
if(~exist(OutputDir,'dir'))
    disp(['Output directory "',OutputDir,'" does not exist, will create it...']);
    mkdir(OutputDir);
end

%% save ClusterEncodeStruct
disp(['Saving "ClusterEncodeStruct.mat" to directory "',OutputDir,'"...']);
save([OutputDir,filesep,'ClusterEncodeStruct.mat'],'ClusterEncodeStruct');

%% output this ClusterEncode Matrix as 4D NIFTI-file
disp(['Writting out "ClusterEncode4D.nii" to directory "',OutputDir,'"...']);
Vout = rmfield(V_ClusteredThreshICs(1),'private'); %make a copy and remove possibly problematic NIFTI field
Vout.fname = [OutputDir,filesep,'ClusterEncode4D.nii'];
if(Vout.dt(1)<16)
    Vout.dt(1) = 16; %just to be save
end
reverseStr = '';
for IndCL = 1:size(ClusterEncode2D,2)
    msg = sprintf(['Writing out Volume/Cluster %0',num2str(ceil(log10(size(ClusterEncode2D,2)))),'d of %0',num2str(ceil(log10(size(ClusterEncode2D,2)))),'d...'], IndCL, size(ClusterEncode2D,2));
    fprintf([reverseStr, msg]); %delete the last message and print current message.
    reverseStr = repmat(sprintf('\b'), 1, length(msg)); %setup the next reverse string
    
    Vout.n(1) = IndCL;
    Data = zeros(size(Mask3D,1),size(Mask3D,2),size(Mask3D,3));
    Data(Mask3D(:)~=0) = ClusterEncode2D(:,IndCL);
    spm_write_vol(Vout,Data);
end

%% Output Uniques summed together, i.e. if there are other values there than 1s (ones) then there are overlaps!
UniqueClusterEncode2D = sum(ClusterEncode2D(:,NClusters==1),2);
disp(['Writting out "UniqueClusterEncode3D.nii" to directory "',OutputDir,'"...']);
Vout = rmfield(V_ClusteredThreshICs(1),'private'); %make a copy and remove possibly problematic NIFTI field
Vout.fname = [OutputDir,filesep,'UniqueClusterEncode3D.nii'];
if(Vout.dt(1)<16)
    Vout.dt(1) = 16; %just to be save
end
Vout.n(1) = 1;
Data = zeros(size(Mask3D,1),size(Mask3D,2),size(Mask3D,3));
Data(Mask3D(:)~=0) = UniqueClusterEncode2D;
spm_write_vol(Vout,Data);

%% Output Overlaps ONLY: first thresholded and then summed together
OverlapsOnlyClusterEncode2D = sum(ClusterEncode2D(:,NClusters>1)~=0,2);
disp(['Writting out "OverlapsOnlyClusterEncode3D.nii" to directory "',OutputDir,'"...']);
Vout = rmfield(V_ClusteredThreshICs(1),'private'); %make a copy and remove possibly problematic NIFTI field
Vout.fname = [OutputDir,filesep,'OverlapsOnlyClusterEncode3D.nii'];
if(Vout.dt(1)<16)
    Vout.dt(1) = 16; %just to be save
end
Vout.n(1) = 1;
Data = zeros(size(Mask3D,1),size(Mask3D,2),size(Mask3D,3));
Data(Mask3D(:)~=0) = OverlapsOnlyClusterEncode2D;
spm_write_vol(Vout,Data);


%% to do
%reanalyse this 4D NIFTI to make unique and overlapping cluster displays

%% Done.
disp(' ');
disp('Done.');
disp(' ');

end

%% subfunction
%% findClustersPerIC
function [ClusterICoverlap, NVoxelClOverlap] = findBestClusterPerIC(OverlapStruct,ICnumClNum,IndCL,OverlapMatThresholds)
% find the best fitting cluster per IC

SelfIC   = ICnumClNum(IndCL,1); %the initial cluster IC should of course be used
SelfClID = ICnumClNum(IndCL,2); %the initial cluster number per IC should of course be used

NICs = max(ICnumClNum(:,1));
ClusterICoverlap = zeros(NICs,1); %init
NVoxelClOverlap  = zeros(NICs,1); %init
for IndIC = 1:NICs
    if(IndIC~=SelfIC)
        %check which one fits best
        AllPossibleClNums      = ICnumClNum(ICnumClNum(:,1)==IndIC,2);
        AllPossibleNVoxOverlap = OverlapStruct.OverlapSize( IndCL,ICnumClNum(:,1)==IndIC);
        AllPossibleMinClOverlap= OverlapStruct.MinClOverlap(IndCL,ICnumClNum(:,1)==IndIC);
        AllPossibleOverlaps    = (AllPossibleNVoxOverlap>=OverlapMatThresholds.NVoxelMin) & (AllPossibleMinClOverlap>=OverlapMatThresholds.MinMinClOverlap);
        if(any(AllPossibleOverlaps~=0))
            AllPossibleClNums = AllPossibleClNums(AllPossibleOverlaps); %use only those that survived the above test.
            [NVoxelClOverlap(IndIC),Ind] = max(AllPossibleNVoxOverlap(AllPossibleOverlaps)); %the possible Number of Voxels for the overlap that survived the above test.
            ClusterICoverlap(IndIC)      = AllPossibleClNums(Ind);
        end
    else
        %assign own
        ClusterICoverlap(IndIC) = SelfClID;
        NVoxelClOverlap(IndIC)  = OverlapStruct.OverlapSize(IndCL,IndCL); %(IndCL,(ICnumClNum(IndCL,1)==SelfIC)&(ICnumClNum(IndCL,2)==SelfClID));
    end
end

end