function [H,OverlapStruct,InfoStruct] = SetupClusterAnalysis(ThreshMapPaths,Thresh,SearchDistance,OutputDir)
% This function sets up the Cluster Analysis of ICA results.
%
% 1.It collects the ICA results, applies a threshold and saves these maps
%   in a 4D NIFTI-file called "ThreshICs.nii", in an output directory that 
%   is created in the directory of the input ICA maps or any given location.
%   This is done by the function "CollectThreshSaveAs4D.m"
% 
% 2.After thresholding, this function will continue by clustering all thresholded maps
%   that have been saved in the previous step.
%   Per Volume the highest significant cluster will be Cluster 1. 
%   Increasing cluster numbers will be associated with decreasing significance value.
%   NB: Statistic values are assumed to be z-, t- or F-statistic values,
%       i.e. larger values are more significant.
%   
% Clustering Scheme:
%   1. Clusters will first be formed via connected voxels here called "superclusters",
%   i.e., all significant voxels that connect with their neighbors form clusters.
%
%   2. These "superclusters" will be further separated by forming clusters based on 
%   the significance values of the contained voxels and a "search distance".
%
%   The "search distance" indicates the distance in mm (MNI-coordinates) that each voxel
%   can "search" through for a local maximum, i.e. another location that is higher in
%   statistical significance. 
%   This is done iteratively until the voxels don't change their association with local maxima any more.
%
%   In other words, a voxel starts out with its local maxima position as its own location,
%   and then looks within the search distance for a location that is higher in significance
%   and sets that as its new local maximum, then the same is repeated from this location
%   until there is no change any more.
%   The default "search distance" is 6 mm, but this does not mean that all clusters will end up having (6mm)^3 in volume!
%   However, local maxima (defining clusters) cannot be closer together than 6 mm or these would converge to one cluster.
%
%   NB: if the "search distance" is set to 0 (zero), then the clusters will be formed by connecting voxels
%       with their neighbors only, i.e. only "superclusters".
%
% 3.Overlap matrices are calculated that can help determining the overlaps that occur
%   and threshold them such that only meaningful overlaps are allowed.
%
%
%
%Usage:
%       [H,OverlapStruct,InfoStruct] = SetupClusterAnalysis(ThreshMapPaths,Thresh,SearchDistance,OutputDir);
%       [H,OverlapStruct,InfoStruct] = SetupClusterAnalysis(ThreshMapPaths,Thresh,SearchDistance); %Determine OutputDir from filepaths in ThreshMapPaths (if possible)  
%       [H,OverlapStruct,InfoStruct] = SetupClusterAnalysis(ThreshMapPaths,Thresh); %set SearchDistance to default value 6 mm AND determine OutputDir from filepaths in ThreshMapPaths (if possible)  
%       [H,OverlapStruct,InfoStruct] = SetupClusterAnalysis(ThreshMapPaths); %Take all values that are non-zero, i.e. Thresh = [eps, -eps]; AND set SearchDistance to default value 6 mm AND determine OutputDir from filepaths in ThreshMapPaths (if possible)  
%       [H,OverlapStruct,InfoStruct] = SetupClusterAnalysis(); %select files manually using spm_select AND take all values that are non-zero, i.e. Thresh = [eps, -eps]; AND set SearchDistance to default value 6 mm AND determine OutputDir from filepaths in ThreshMapPaths (if possible)  
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (01.11.2016): initial implementation based on previous tests

%% Check inputs
%ThreshMapPaths
if(~exist('ThreshMapPaths','var'))
    ThreshMapPaths = cellstr(spm_select(inf,'image','Select IC maps NIFTI-files...'));
else
    if(isempty(ThreshMapPaths))
        ThreshMapPaths = cellstr(spm_select(inf,'image','Select IC maps NIFTI-files...'));
    else
        if(~iscellstr(ThreshMapPaths))
            if(ischar(ThreshMapPaths))
                ThreshMapPaths = cellstr(ThreshMapPaths);
            else
                error('ThreshMapPaths has to be a cellstring or convertible char-array, indicating the paths of the NIFTI files!');
            end
        end
    end
end

%Thresh
ThreshIsEmpty = 1; %init true (this will be used to keep track of the initial state of thresholds
if(~exist('Thresh','var'))
    Thresh = [eps, -eps];
else
    if(isempty(Thresh))
        Thresh = [eps, -eps];
    else
        ThreshIsEmpty = 0;
        if((sign(Thresh(1))==-1) && ((sign(Thresh(2))==+1))) %need to flip it
            Thresh = Thresh(2:-1:1);
        else
            if(((sign(Thresh(1))==-1) && (sign(Thresh(2))==-1)) || ((sign(Thresh(1))==+1) && (sign(Thresh(2))==+1)))
                error('Thresh should have a positive and a negative value or empty if none used.');
            end                
        end
    end
end
disp(['Will use Thresholds [',num2str(Thresh(:)'),'] on IC maps.']);

%SearchDistance NB: no messages here because Cluster4DStatsNIFTI.m outputs messages.
if(~exist('SearchDistance','var'))
    %disp('No SearchDistance entered, will use default search distance of 6 mm.');
    SearchDistance = 6; %mm
else
    if(isempty(SearchDistance))
        %disp('Will cluster stats NIFTIs only into "superclusters", i.e., via connected voxels and will NOT separate these further (SearchDistance==0)!');
        SearchDistance = 0;
    end
end
if(SearchDistance<0)
    error('SearchDistance must be >=0! (==0 mean it is not used, only "superclusters" are formed)');    
end

%OutputDir
if(~exist('OutputDir','var'))
    OutputDir = unique(cellfun(@fileparts,ThreshMapPaths,'UniformOutput',false));
    if(length(OutputDir)~=1)
        error('OutputDir automatic generation does not work if multiple paths are possible (as extracted from the input ThreshMapPaths)! Please specify the output directory in that case.');
    else
        if(ThreshIsEmpty)
            OutputDir = [OutputDir{1},filesep,'ClusterAnalysis'];
        else
            OutputDir = [OutputDir{1},filesep,'ClusterAnalysis_',num2str(Thresh(1)),'_',num2str(Thresh(2))];
        end
    end
else
    if(isempty(OutputDir))
        OutputDir = unique(cellfun(@fileparts,ThreshMapPaths,'UniformOutput',false));
        if(length(OutputDir)~=1)
            error('OutputDir automatic generation does not work if multiple paths are possible (as extracted from the input ThreshMapPaths)! Please specify the output directory in that case.');
        else
            if(ThreshIsEmpty)
                OutputDir = [OutputDir{1},filesep,'ClusterAnalysis'];
            else
                OutputDir = [OutputDir{1},filesep,'ClusterAnalysis_',num2str(Thresh(1)),'_',num2str(Thresh(2))];
            end
        end
    end
end

%% create ThreshICs 4D-NIFTI in OutputDir
OutPathThreshICs = [OutputDir,filesep,'ThreshICs.nii'];
InfoStruct.CollectionThreshStruct = CollectThreshSaveAs4D(ThreshMapPaths,Thresh,OutPathThreshICs);

%% cluster ThreshICs and write 4D-NIFTI file
OutPath4DClusterNIFTI = [OutputDir,filesep,'Clustered',num2str(SearchDistance),'mmThreshICs.nii'];
InfoStruct.ClusteringStruct = Cluster4DStatsNIFTI(OutPathThreshICs,SearchDistance,OutPath4DClusterNIFTI);

%% save InfoStruct
disp(['Saving InfoStruct.mat to directory "',OutputDir,'".']);
save([OutputDir,filesep,'InfoStruct.mat'],'InfoStruct');

%% GenerateOverlapMatrices
OverlapStruct = GenerateOverlapMatrices(OutputDir);

%% DisplayOverlapStructData
H = DisplayOverlapStructData(OverlapStruct);

%% Done.
disp(' ');
disp('Done.');
disp(' ');

end

