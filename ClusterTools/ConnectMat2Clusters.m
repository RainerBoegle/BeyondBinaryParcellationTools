function [ClusterNr_vox, ClusterSizes] = ConnectMat2Clusters(ConnectMat,StatsVals_vox)
% This function can be used to turn a connectedness matrix into a clustering of voxels simply via
% their connectedness.
%
% WARNING: 
%          This is usually a bad idea, as connectedness does not separate clusters well.
%          You better know what you are doing!
% 
%
%Usage:
%      [ClusterNr_vox, ClusterSizes] = ConnectMat2Clusters(ConnectMat,StatsVals_vox); %Pick highest StatsVal-Voxel and all connected to it and assign them as Cluster 1,
%                                                                                     %for the remaining look again for highest and connected and assign them to Cluster 2, 
%                                                                                     %continue until no voxels left.
%      [ClusterNr_vox, ClusterSizes] = ConnectMat2Clusters(ConnectMat); %cluster number 1 is for the first voxel on the list and all connected to it, 
%                                                                       %then number 2 is the next on the list that is not connected to the first cluster and so on until no voxels left.
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment(07.November.2015): initial implementation based on test script.

%% Check inputs
if(isempty(ConnectMat))
    error('ConnectMat can not be empty!');
else
    NVoxels = size(ConnectMat,1);
end

try
    if(isempty(StatsVals_vox))
        StatsVals_vox = (NVoxels:-1:1)'; %This way the first voxel is "the most significant one". 
        disp('Can not use "StatsVals_vox". Clusters will start with first voxels and connected and continue from there with the remaining, till no voxels left.');
    else
        if(length(StatsVals_vox)~=NVoxels)
            error('"StatsVals_vox" need to contain as many values(==Nvoxels) as (corresponding) number of voxels in connectedness matrix!');
        end
    end
catch CATCH_StatsVals_vox
    StatsVals_vox = (NVoxels:-1:1)'; 
    disp('NB: No StatsVals --> Clusters will start with first voxels and connected and continue from there with the remaining, till no voxels left.');
end

%% pick a Voxel (based on StatsVals_vox, -whatever that is at this point, either real stats-vals or just a list "NVoxels:-1:1", see above.)
%% then get all that are connected to this one
%% assign current cluster number to these voxels.
%% remove from list and repeat until done.
disp('Assigning clusters on basis of connectedness matrix...');
ClusterNr_vox       = zeros(NVoxels,1); %init as zeros, i.e. empty.
CurrClusterNum      = 1; %counter
CurrClusterSize     = []; %init empty
[Stmp,VoxelIdxToDo] = sort(StatsVals_vox,1,'descend'); clear Stmp %from "most significant" to "least significant"
while(~isempty(VoxelIdxToDo))
    CurrStartIdx = VoxelIdxToDo(1); %pick a voxel to start from.
    ClusterNr_vox(ConnectMat(CurrStartIdx,:)~=0) = CurrClusterNum; %assign current cluster number to voxels that are connected to the one that was just picked.
    IdxToRemove = find(ConnectMat(CurrStartIdx,:)~=0);
    for IndRemove = 1:length(IdxToRemove)
        VoxelIdxToDo(VoxelIdxToDo==IdxToRemove(IndRemove)) = []; %remove those that are connected and have just been found from list.
    end
    CurrClusterSize = [CurrClusterSize; length(IdxToRemove)];
    CurrClusterNum = CurrClusterNum + 1; %increase cluster number for next iteration.
end

%% do some checks
if(max(ClusterNr_vox)>1)
    disp([num2str(max(ClusterNr_vox)),' Clusters have been assigned. (ClusterSizes are [',num2str(CurrClusterSize'),'] Voxels for ',num2str(NVoxels),' Voxels in total.)']);
else
    if(all(ClusterNr_vox==0))
        error('No clusters have been assigned!');
    else
        disp(['WARNING: only ONE cluster has been assigned for ',num2str(NVoxels),' Voxels.']);
    end
end

%% get cluster sizes
ClusterSizes = zeros(max(ClusterNr_vox),1);
for IndCl = 1:length(ClusterSizes)
    ClusterSizes(IndCl) = length(find(ClusterNr_vox==IndCl));
end

%% done.
disp('...done.');
disp(' ');

end

