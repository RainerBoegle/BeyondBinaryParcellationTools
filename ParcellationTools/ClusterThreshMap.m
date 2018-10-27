function [ClusterMapPath,VolCL,LocMaxStruct,XYZmm,StatsVals,XYZvox,VolStats] = ClusterThreshMap(ThreshMapPath,Thresh,UseConnectMat,SearchDist_mm,MinClusterSize)
% This function can be used to get thresholded data from a stats NIFTI-file and cluster the thresholded stats NIFTI-file.
% Thresh indicates positive and negative threshold.
% UseConnectMat is 0or1 indicating if clusters should be initially formed via connected voxels.
% SearchDist_mm is the distance that each voxel can "move" (if UseConnectMat==1 then they can only "move" within the connected cluster)
% to the next "higher"/more significant voxel that is in distance SearchDist_mm. 
% This is then repeated until convergence.
% NB: each cluster is tested for its size using "MinClusterSize" and clusters smaller than "MinClusterSize" will be thrown out.
%
%Usage:
%      [ClusterMapPath,VolCL,LocMaxStruct,XYZmm,StatsVals,XYZvox,VolStats] = ClusterThreshMap(ThreshMapPath,Thresh,UseConnectMat,SearchDist_mm,MinClusterSize);
%      [ClusterMapPath,VolCL,LocMaxStruct,XYZmm,StatsVals,XYZvox,VolStats] = ClusterThreshMap(ThreshMapPath,[2.3 -Inf],1,8,0); %load map and threshold it with Stats>2.3 & Stats<-Inf, i.e. don't get negative stats; initially create clusters as connected regions (UseConnectMat) and then find local maxima by iteratively moving voxels to "higherst" voxel in a distance of 8mm (SearchDist_mm) until convergence.
%      [ClusterMapPath,VolCL,LocMaxStruct,XYZmm,StatsVals,XYZvox,VolStats] = ClusterThreshMap(ThreshMapPath,[0 0],     1,8,0); %the same as above but assume map is already thresholded.
%      [ClusterMapPath,VolCL,LocMaxStruct,XYZmm,StatsVals,XYZvox,VolStats] = ClusterThreshMap(ThreshMapPath,[0 0],     1,8,6); %the same as above but remove clusters smaller than 6 voxels.
%
%V1.5
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.5: (22.03.2016): add cluster size exclusion. V1.0: (16.03.2016): initial implementation based on previous tests


%% Check inputs
%MinClusterSize
if(~exist('MinClusterSize','var'))
    MinClusterSize = 0; %default: DO NOT EXCLUDE CLUSTERS based on size.
else
    if(isempty(MinClusterSize))
        MinClusterSize = 0; %default: DO NOT EXCLUDE CLUSTERS based on size.
    else
        if(MinClusterSize<0)
            error(['A minimum cluster size of ',num2str(MinClusterSize),' is not valid!']);
        end
    end
end
if(MinClusterSize~=0)
    disp(['Will exclude clusters of size smaller than ',num2str(MinClusterSize),'...']);
else
    disp('Warning! Will not exclude clusters based on their size. This may lead to very small clusters!');
end

%% load inputs
disp('Loading input map');
[XYZmm,StatsVals,XYZvox,VolStats,VolInd] = GetDataFromMap(ThreshMapPath,Thresh);
if(~isempty(VolInd))
    VolStr = ['_',num2str(VolInd,['%0',num2str(max([ceil(log10(VolInd)); 2])),'d'])];
else 
    VolStr = '';
end

%% Do clustering
disp('Clustering thresholded map');
if(UseConnectMat)
    ConnectMat   = cell(2,1);
    LocMaxStruct = cell(2,1);
    ConnMatStr   = ' & restricting to connected voxels';
    %positive stats vals
    if(~isempty(XYZvox{1}))
        disp('Positive clusters will be examined');
        disp('Determining connectedness...');
        ConnectMat{1} = FindConnectedVoxels(XYZvox{1});
        if(MinClusterSize~=0) %exclude first time
            [ClusterNr_vox, ClusterSizes] = ConnectMat2Clusters(ConnectMat{1});
            IndClExclude = find(ClusterSizes<MinClusterSize);
            while(~isempty(IndClExclude))
                ClusterNr_vox_Exclude = ClusterNr_vox;
                disp([num2str(length(IndClExclude)),' Clusters will be removed due to cluster size < ',num2str(MinClusterSize),'. (after connected clusters estimation)']);
                for IndClEx = 1:length(IndClExclude)
                    CurrClInd = IndClExclude(IndClEx);
                    ClusterNr_vox_Exclude(ClusterNr_vox_Exclude==CurrClInd) = 0; %mark for removal
                end
                XYZvox{1}   = XYZvox{1}(ClusterNr_vox_Exclude~=0,:); %reassign
                StatsVals{1}= StatsVals{1}(ClusterNr_vox_Exclude~=0);
                XYZmm{1}    = XYZmm{1}(ClusterNr_vox_Exclude~=0,:); %reassign
                
                disp('Re-Determining connectedness...');
                ConnectMat{1}= FindConnectedVoxels(XYZvox{1});
                [ClusterNr_vox, ClusterSizes] = ConnectMat2Clusters(ConnectMat{1});
                IndClExclude = find(ClusterSizes<MinClusterSize);
            end
        end
        disp('Reclustering connected voxels...');
        LocMaxStruct{1} = FindAllLocMax(XYZmm{1},StatsVals{1},SearchDist_mm,{'VoxelConnections',ConnectMat{1}});
        if(MinClusterSize~=0) %exclude second time
            IndClExclude = find(LocMaxStruct{1}.LocMaxClusterSize<MinClusterSize);
            while(~isempty(IndClExclude))
                ClusterNr_vox_Exclude = LocMaxStruct{1}.ClusterNo;
                disp([num2str(length(IndClExclude)),' Clusters will be removed due to cluster size < ',num2str(MinClusterSize),'. (after FindAllLocMax estimation)']);
                for IndClEx = 1:length(IndClExclude)
                    CurrClInd = IndClExclude(IndClEx);
                    ClusterNr_vox_Exclude(ClusterNr_vox_Exclude==CurrClInd) = 0; %mark for removal
                end
                XYZvox{1}   = XYZvox{1}(ClusterNr_vox_Exclude~=0,:); %reassign
                StatsVals{1}= StatsVals{1}(ClusterNr_vox_Exclude~=0);
                XYZmm{1}    = XYZmm{1}(ClusterNr_vox_Exclude~=0,:); %reassign
                
                disp('Re-Determining connectedness...');
                ConnectMat{1}= FindConnectedVoxels(XYZvox{1});
                
                disp('Re-Reclustering connected voxels...');
                LocMaxStruct{1} = FindAllLocMax(XYZmm{1},StatsVals{1},SearchDist_mm,{'VoxelConnections',ConnectMat{1}});
                IndClExclude = find(LocMaxStruct{1}.LocMaxClusterSize<MinClusterSize);
            end
        end
        disp(' ');
    else
        ConnectMat{1}   = [];
        LocMaxStruct{1} = [];
    end
    %negative stats vals
    if(~isempty(XYZvox{2}))
        disp('Negative clusters will be examined');
        disp('Determining connectedness...');
        ConnectMat{2} = FindConnectedVoxels(XYZvox{2});
        if(MinClusterSize~=0) %exclude first time
            [ClusterNr_vox, ClusterSizes] = ConnectMat2Clusters(ConnectMat{2});
            IndClExclude = find(ClusterSizes<MinClusterSize);
            while(~isempty(IndClExclude))
                ClusterNr_vox_Exclude = ClusterNr_vox;
                disp([num2str(length(IndClExclude)),' Clusters will be removed due to cluster size < ',num2str(MinClusterSize),'. (after connected clusters estimation)']);
                for IndClEx = 1:length(IndClExclude)
                    CurrClInd = IndClExclude(IndClEx);
                    ClusterNr_vox_Exclude(ClusterNr_vox_Exclude==CurrClInd) = 0; %mark for removal
                end
                XYZvox{2}   = XYZvox{2}(ClusterNr_vox_Exclude~=0,:); %reassign
                StatsVals{2}= StatsVals{2}(ClusterNr_vox_Exclude~=0);
                XYZmm{2}    = XYZmm{2}(ClusterNr_vox_Exclude~=0,:); %reassign
                
                disp('Re-Determining connectedness...');
                ConnectMat{2}= FindConnectedVoxels(XYZvox{2});
                [ClusterNr_vox, ClusterSizes] = ConnectMat2Clusters(ConnectMat{2});
                IndClExclude = find(ClusterSizes<MinClusterSize);
            end
        end
        disp('Reclustering connected voxels...');
        LocMaxStruct{2} = FindAllLocMax(XYZmm{2},StatsVals{2},SearchDist_mm,{'VoxelConnections',ConnectMat{2}});
        if(MinClusterSize~=0) %exclude second time
            IndClExclude = find(LocMaxStruct{2}.LocMaxClusterSize<MinClusterSize);
            while(~isempty(IndClExclude))
                ClusterNr_vox_Exclude = LocMaxStruct{2}.ClusterNo;
                disp([num2str(length(IndClExclude)),' Clusters will be removed due to cluster size < ',num2str(MinClusterSize),'. (after FindAllLocMax estimation)']);
                for IndClEx = 1:length(IndClExclude)
                    CurrClInd = IndClExclude(IndClEx);
                    ClusterNr_vox_Exclude(ClusterNr_vox_Exclude==CurrClInd) = 0; %mark for removal
                end
                XYZvox{2}   = XYZvox{2}(ClusterNr_vox_Exclude~=0,:); %reassign
                StatsVals{2}= StatsVals{2}(ClusterNr_vox_Exclude~=0);
                XYZmm{2}    = XYZmm{2}(ClusterNr_vox_Exclude~=0,:); %reassign
                
                disp('Re-Determining connectedness...');
                ConnectMat{2}= FindConnectedVoxels(XYZvox{2});
                
                disp('Re-Reclustering connected voxels...');
                LocMaxStruct{2} = FindAllLocMax(XYZmm{2},StatsVals{2},SearchDist_mm,{'VoxelConnections',ConnectMat{2}});
                IndClExclude = find(LocMaxStruct{2}.LocMaxClusterSize<MinClusterSize);
            end
        end
        disp(' ');
    else
        ConnectMat{2}   = [];
        LocMaxStruct{2} = [];
    end
else
    LocMaxStruct = cell(2,1);
    ConnMatStr   = '';
    %positive stats vals
    if(~isempty(XYZvox{1}))
        disp('Positive clusters will be examined');
        disp('Clustering voxels...');
        LocMaxStruct{1} = FindAllLocMax(XYZmm{1},StatsVals{1},SearchDist_mm);
        if(MinClusterSize~=0) %exclude first time
            IndClExclude = find(LocMaxStruct{1}.LocMaxClusterSize<MinClusterSize);
            while(~isempty(IndClExclude))
                ClusterNr_vox_Exclude = LocMaxStruct{1}.ClusterNo;
                disp([num2str(length(IndClExclude)),' Clusters will be removed due to cluster size < ',num2str(MinClusterSize),'. (after FindAllLocMax estimation)']);
                for IndClEx = 1:length(IndClExclude)
                    CurrClInd = IndClExclude(IndClEx);
                    ClusterNr_vox_Exclude(ClusterNr_vox_Exclude==CurrClInd) = 0; %mark for removal
                end
                XYZvox{1}   = XYZvox{1}(ClusterNr_vox_Exclude~=0,:); %reassign
                StatsVals{1}= StatsVals{1}(ClusterNr_vox_Exclude~=0);
                XYZmm{1}    = XYZmm{1}(ClusterNr_vox_Exclude~=0,:); %reassign
                
                disp('Re-Reclustering connected voxels...');
                LocMaxStruct{1} = FindAllLocMax(XYZmm{1},StatsVals{1},SearchDist_mm);
                IndClExclude = find(LocMaxStruct{1}.LocMaxClusterSize<MinClusterSize);
            end
        end
        disp(' ');
    else
        LocMaxStruct{1} = [];
    end
    %negative stats vals
    if(~isempty(XYZvox{2}))
        disp('Negative clusters will be examined');
        disp('Clustering voxels...');
        LocMaxStruct{2} = FindAllLocMax(XYZmm{2},StatsVals{2},SearchDist_mm);
        if(MinClusterSize~=0) %exclude first time
            IndClExclude = find(LocMaxStruct{2}.LocMaxClusterSize<MinClusterSize);
            while(~isempty(IndClExclude))
                ClusterNr_vox_Exclude = LocMaxStruct{2}.ClusterNo;
                disp([num2str(length(IndClExclude)),' Clusters will be removed due to cluster size < ',num2str(MinClusterSize),'. (after FindAllLocMax estimation)']);
                for IndClEx = 1:length(IndClExclude)
                    CurrClInd = IndClExclude(IndClEx);
                    ClusterNr_vox_Exclude(ClusterNr_vox_Exclude==CurrClInd) = 0; %mark for removal
                end
                XYZvox{2}   = XYZvox{2}(ClusterNr_vox_Exclude~=0,:); %reassign
                StatsVals{2}= StatsVals{2}(ClusterNr_vox_Exclude~=0);
                XYZmm{2}    = XYZmm{2}(ClusterNr_vox_Exclude~=0,:); %reassign
                
                disp('Re-Reclustering connected voxels...');
                LocMaxStruct{2} = FindAllLocMax(XYZmm{2},StatsVals{2},SearchDist_mm);
                IndClExclude = find(LocMaxStruct{2}.LocMaxClusterSize<MinClusterSize);
            end
        end
        disp(' ');
    else
        LocMaxStruct{2} = [];
    end
end

%% write out cluster maps
[BaseDir,FNameOrg] = fileparts(ThreshMapPath);

if(UseConnectMat)
    ConnStr = 'Conn';
else
    ConnStr = '';
end
if(MinClusterSize~=0)
    MinClusterSizeStr = ['MinCl',num2str(MinClusterSize)];
else
    MinClusterSizeStr = '';
end
ClusterMapPath= cell(2,1);
VolCL         = cell(2,1);
for Ind = 1:2
    if(isstruct(LocMaxStruct{Ind})&&~isempty(LocMaxStruct{Ind}))
        if(Ind==1)
            disp('Writing out POSITIVE STATS CLUSTERS BEFORE rethresholding.');
            StatsStr = 'PosStats';
        else
            disp('Writing out NEGATIVE STATS CLUSTERS BEFORE rethresholding.');
            StatsStr = 'NegStats';
        end
        
        %Clusters are separated by stats sign (Pos/Neg)
        OutputPath          = [BaseDir,filesep,'C',num2str(SearchDist_mm),'mm',ConnStr,MinClusterSizeStr,StatsStr,'_',FNameOrg,VolStr,'.nii'];
        [ClusterMapPath{Ind},VolCL{Ind}] = LocMaxClusters2NIFTI(LocMaxStruct{Ind}.ClusterNo,XYZvox{Ind},VolStats,OutputPath);
    end
end

%% done
disp(['Done clustering (search dist = ',num2str(SearchDist_mm),'mm',ConnMatStr,').']);
disp('DONE.');
disp(' ');

end