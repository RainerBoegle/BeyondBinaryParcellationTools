function ClusteringStruct = Cluster4DStatsNIFTI(Thresh4DStatsNIFTIPath,SearchDistance,OutPath4DClusterNIFTI)
% This function clusters all volumes of a 4D-NIFTI file assuming that it is a statistic file
% containing z-, t- or F-statistic values, i.e., larger means more significant. 
% Negative and positive values are treated separately, only the amplitude matters. 
%
% NB: It is assumed that the volumes are already thresholded, 
%     i.e. only non-zero values are statistically significant.
%
%Clustering Scheme:
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
%
%
%Usage:
%       ClusteringStruct = Cluster4DStatsNIFTI(Thresh4DStatsNIFTIPath,SearchDistance,OutPath4DClusterNIFTI);
%       ClusteringStruct = Cluster4DStatsNIFTI(Thresh4DStatsNIFTIPath,SearchDistance); %Determine OutPath4DClusterNIFTI from filepaths in Thresh4DStatsNIFTIPath (if possible).  
%       ClusteringStruct = Cluster4DStatsNIFTI(Thresh4DStatsNIFTIPath); %set SearchDistance to default value 6 mm AND determine OutPath4DClusterNIFTI from filepaths in Thresh4DStatsNIFTIPath (if possible).  
%       ClusteringStruct = Cluster4DStatsNIFTI(); %select 4D-NIFTI manually, set SearchDistance to default value 6 mm AND determine OutPath4DClusterNIFTI from filepaths in Thresh4DStatsNIFTIPath as selected by user (if possible).  
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (01.11.2016): initial implementation based on previous tests

%% Check inputs
if(~exist('Thresh4DStatsNIFTIPath','var'))
    Thresh4DStatsNIFTIPath = cellstr(spm_select(inf,'image','Select 4D-NIFTI file with all volumes...'));
else
    if(isempty(Thresh4DStatsNIFTIPath))
        Thresh4DStatsNIFTIPath = cellstr(spm_select(inf,'image','Select 4D-NIFTI file with all volumes...'));
    else
        if(~iscellstr(Thresh4DStatsNIFTIPath))
            if(ischar(Thresh4DStatsNIFTIPath))
                Thresh4DStatsNIFTIPath = cellstr(Thresh4DStatsNIFTIPath);
            else
                error('Thresh4DStatsNIFTIPath has to be a string pointing to the Statistic NIFTI file that should be clustered.');
            end
        end
    end
end

%SearchDistance
if(~exist('SearchDistance','var'))
    disp('No SearchDistance entered, will use default search distance of 6 mm.');
    SearchDistance = 6; %mm
else
    if(isempty(SearchDistance))
        disp('Will cluster stats NIFTIs only into "superclusters", i.e., via connected voxels and will NOT separate these further (SearchDistance==0)!');
        SearchDistance = 0;
    end
end
if(SearchDistance==0)
    disp('Will cluster stats NIFTIs only into "superclusters", i.e., via connected voxels, and will NOT separate these further (SearchDistance==0)!');
else
    if(SearchDistance<0)
        error('SearchDistance must be >=0! (==0 mean it is not used, only "superclusters" are formed)');
    else
        disp( 'Will cluster by creating "superclusters", i.e. connected voxels,');
        disp( 'and then try to separate these into smaller clusters,');
        disp(['using a search distance of ',num2str(SearchDistance),' mm.']);
    end
end 

%OutPath4DClusterNIFTI
if(~exist('OutPath4DClusterNIFTI','var')) %need to produce it via other paths
    [OutPath4DClusterNIFTI,OutFNames] = unique(cellfun(@fileparts,Thresh4DStatsNIFTIPath,'UniformOutput',false));
    if(length(OutPath4DClusterNIFTI)~=1)
        error('OutPath4DClusterNIFTI automatic generation does not work if multiple paths are possible (as extracted from the input Thresh4DStatsNIFTIPath)! Please specify the output path in that case.');
    else
        OutPath4DClusterNIFTI = [OutPath4DClusterNIFTI{1},filesep,'Clustered',num2str(SearchDistance),'mm',OutFNames{1},'.nii'];
        clear OutFNames
    end
else
    if(isempty(OutPath4DClusterNIFTI))
        [OutPath4DClusterNIFTI,OutFNames] = unique(cellfun(@fileparts,Thresh4DStatsNIFTIPath,'UniformOutput',false));
        if(length(OutPath4DClusterNIFTI)~=1)
            error('OutPath4DNIFTI automatic generation does not work if multiple paths are possible (as extracted from the input Thresh4DStatsNIFTIPath)! Please specify the output path in that case.');
        else
            OutPath4DClusterNIFTI = [OutPath4DClusterNIFTI{1},filesep,'Clustered',num2str(SearchDistance),'mm',OutFNames{1},'.nii'];
            clear OutFNames
        end
    else
        if(~ischar(OutPath4DClusterNIFTI))
            if(~iscellstr(OutPath4DClusterNIFTI))
                error('OutPath4DNIFTI must be a string indicating the output path for the 4D-NIFTI file!');
            else
                OutPath4DClusterNIFTI = OutPath4DClusterNIFTI{1};
            end
        end
    end
end

%% get the data
Vols = spm_vol(Thresh4DStatsNIFTIPath);
[StatsVals,XYZvox,XYZmm] = GetDataFromNII(Vols);

%% do clustering
ClNumsVox  = cell(length(StatsVals),1);
for IndVol = 1:length(StatsVals)
    disp(['Clustering Volume ',num2str(IndVol,['%0',num2str(ceil(log10(length(StatsVals)))),'d']),' of ',num2str(length(StatsVals),['%0',num2str(ceil(log10(length(StatsVals)))),'d']),'...']);
    
    ClNumsVox{IndVol} = DoClustering(StatsVals{IndVol},XYZvox{IndVol},XYZmm{IndVol},SearchDistance);
    disp(' ');
end

%% create output volume 
if(iscell(Vols))
    if(isstruct(Vols{1}))
        Vout = rmfield(Vols{1},'private'); %prepare by removing any NIFTI attributes that could cause trouble.
    end
else
    if(isstruct(Vols))
        Vout = rmfield(Vols,'private'); %prepare by removing any NIFTI attributes that could cause trouble.
    end
end
if(length(Vout)>1)
    Vout = Vout(1);
end
if(Vout.dt(1)<16)
    Vout.dt(1) = 16; %make sure this is encoded at least as a 32bit float.
end
Vout.fname = OutPath4DClusterNIFTI;

%% check if output dir exists otherwise create it
[OutDir,OutFName,OutExt] = fileparts(OutPath4DClusterNIFTI);
if(~exist(OutDir,'dir'))
    disp(['Directory "',OutDir,'" does not exist, will create it...']);
    mkdir(OutDir);
end

%% write out data
disp(['Writing out "',OutFName,OutExt,'" to directory "',OutDir,'"...']);
reverseStr = ''; %init
for IndVol = 1:length(ClNumsVox)
    Y = zeros(Vout.dim(1),Vout.dim(2),Vout.dim(3));
    CurrClNumsVox = ClNumsVox{IndVol};
    CurrXYZvox= XYZvox{IndVol};
    for IndVox = 1:length(CurrClNumsVox)
        Y(CurrXYZvox(IndVox,1),CurrXYZvox(IndVox,2),CurrXYZvox(IndVox,3)) = CurrClNumsVox(IndVox);
    end
    Vout.n(1) = IndVol;
    spm_write_vol(Vout,Y);
    
    msg = sprintf(['Writing out Volume %0',num2str(ceil(log10(length(ClNumsVox)))),'d of %0',num2str(ceil(log10(length(ClNumsVox)))),'d...'], IndVol, length(ClNumsVox));
    fprintf([reverseStr, msg]); %delete the last message and print current message.
    reverseStr = repmat(sprintf('\b'), 1, length(msg)); %setup the next reverse string
end
msg = sprintf(['Writing out Volume %0',num2str(ceil(log10(length(ClNumsVox)))),'d of %0',num2str(ceil(log10(length(ClNumsVox)))),'d done.'], IndVol, length(ClNumsVox));
fprintf([reverseStr, msg]); %delete the last message and print current message.
disp(' ');

%% save ClusteringStruct to OutputDir
ClusteringStruct = struct('Thresh4DStatsNIFTIPath',{Thresh4DStatsNIFTIPath},'SearchDistance',SearchDistance,'OutPath4DClusterNIFTI',OutPath4DClusterNIFTI);
disp(['Saveing ClusteringStruct.mat to directory "',OutDir,'".']);
save([OutDir,filesep,'ClusteringStruct.mat'],'ClusteringStruct');

%% Done.
disp(' ');
disp('Done.');
disp(' ');

end


%% subfunction
function [ClusterNumbersPerVoxel] = DoClustering(StatsVals,XYZvox,XYZmm,SearchDist_mm)
% This function does the clustering and returns the ClusterNumbersPerVoxel

ClusterNumbersPerVoxel = zeros(length(StatsVals),1); %init
if(~isempty(StatsVals(StatsVals>0)))
    disp('Positive clusters will be examined');
    disp('Determining "superclusters", i.e. voxel connectedness...');
    ConnectMat   = FindConnectedVoxels(XYZvox(StatsVals>0,:));
    if(SearchDist_mm>0)
        disp('Reclustering "superclusters" (connected voxels)...');
        LocMaxStruct = FindAllLocMax(XYZmm(StatsVals>0,:),StatsVals(StatsVals>0),SearchDist_mm,{'VoxelConnections',ConnectMat});
        ClusterNumbersPerVoxel(StatsVals>0) = LocMaxStruct.ClusterNo; %(NVoxelx1) The (final) cluster number per input voxel.
        clear ConnectMat LocMaxStruct
    else
        ClusterNumbersPerVoxel(StatsVals>0) = ConnectMat2Clusters(ConnectMat,StatsVals(StatsVals>0));
        clear ConnectMat
    end
    disp('done.');
    disp(' ');
end

Offset = max(ClusterNumbersPerVoxel(:));
%negative stats vals
if(~isempty(StatsVals(StatsVals<0)))
    disp('Negative clusters will be examined');
    disp('Determining "superclusters", i.e. voxel connectedness...');
    ConnectMat   = FindConnectedVoxels(XYZvox(StatsVals<0,:));
    if(SearchDist_mm>0)
        disp('Reclustering "superclusters" (connected voxels)...');
        LocMaxStruct = FindAllLocMax(XYZmm(StatsVals<0,:),StatsVals(StatsVals<0),SearchDist_mm,{'VoxelConnections',ConnectMat});
        ClusterNumbersPerVoxel(StatsVals<0) = Offset+LocMaxStruct.ClusterNo; %(NVoxelx1) The (final) cluster number per input voxel.
        clear ConnectMat LocMaxStruct
    else
        ClusterNumbersPerVoxel(StatsVals<0) = Offset+ConnectMat2Clusters(ConnectMat,StatsVals(StatsVals<0));
        clear ConnectMat
    end
    disp('done.');
    disp(' ');
end

end