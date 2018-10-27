function [ReThreshMapPath,ClusterMapPath,StatsValsFinal,XYZmmFinal,XYZvoxFinal,StatsVals,XYZmm,XYZvox,StatsValsORG,XYZmmORG,XYZvoxORG] = ClusterReThresh(ThreshMapPath,SearchDist_mm,UseConnectMat,RetainPercentage,RetainMode,varargin)
% This function can be used to rethreshold functional maps based on clustering of the maps for local
% maxima and the thresholding of each cluster to retain a certain percentage of voxels close to respective
% maximum value.
% In other words, the map can be split into clusters around local maxima, the voxels of these local
% maxima are then analysed for their respective maximum and only the voxels that are in the top X%
% near that maximum value (of the respective cluster) are kept, the others are dropped.
%
% Inputs:
%         ThreshMapPath     (string)  --> Path to the NIFTI-file containing the functional map
%                                        (should be thresholded already, i.e. the nonzero voxels should be the significant ones) 
%                                        (In case you want this function to threshold your map
%         SearchDist_mm     (double)  --> Distance in mm for clustering: it is assumed that the the clusters have AT LEAST this distance. 
%                                        (Also called the search distance, because each voxels "searches" for a local maximum
%                                         to "belong to" in maximally this distance, which is iterated until convergence is reached.
%                                         This also means that some voxels may find a maximum further than this distance away, but that
%                                         indicates that there is a steady gradient for the voxels to follow, i.e. no local maxima inbetween.)
%         UseConnectMat     (0 or 1)  --> Restrict "movement/search" of voxels for a local maxima in the "search distance to voxels
%                                         that are reachable via connected neighbors. This means that clusters can only form from connected voxels.
%         RetainPercentage  (0<P<=1)  --> Indicates the top P% percent voxels that should be kept, 
%                                         See "RetainMode".
%         RetainMode        (string)  --> Either 'Thresh' (or 'T' or 'Threshold') OR 'Voxels' (or 'V' or 'Vox')
%                                         'Thresh': percentage relative to the maximum of the cluster local maxima and then taking the voxels above this value.
%                                         'Voxels': percentage of VOXELS in the cluster of the local maxima, i.e. from local maximum going downward until "RetainPercentage" voxels (sorted descending).
%
% Xtra-Inputs:
%      {'thresh',[PosLimit NegLimit]} --> thresholding of initial map BEFORE CLUSTERING! {'thresh',[4,-3]} to threshold the initial map (before clustering) with thresholds >4 & <-3.
%      {'CLsize',ClusterSizeThresh}   --> Threshold for number of voxels. Throwing out all that are <=ClusterSizeThresh.
%
% Outputs:
%         ReThreshMapPath   (string)  --> The rethresholded map. Same directory as for Input ThreshMapPath; Filename + prefix(ReThresh"RetainPercentage*100"C"SearchDist"mm""Conn"_)
%                                        (e.g. Filename = 'ABCD' (Cluster 8mm + ConnectMat + 0.3 RetainPercentage) --> OutName = 'ReThresh30C8mmConn_ABCD') 
%         ClusterMapPath    (string)  --> The clusters used for rethresholding. Same directory as for Input ThreshMapPath; Filename + prefix("C"SearchDist"mm""Conn"_)
%                                        (e.g. Filename = 'ABCD' (Cluster 8mm + ConnectMat) --> OutName = 'C8mmConnPos/NegStats_ABCD')
%
%
%Usage:
%        [ReThreshMapPath,ClusterMapPath] = ClusterReThresh(ThreshMapPath,SearchDist_mm,UseConnectMat,RetainPercentage,{'comm-string',values});
%        [ReThreshMapPath,ClusterMapPath] = ClusterReThresh(ThreshMapPath,8,1,0.3); %take (already) thresholded map, cluster it with search distance of 8mm and allow cluster formation only for connected voxels; -then take the top 30% (0.3) of voxels for each cluster.
%        [ReThreshMapPath,ClusterMapPath] = ClusterReThresh(ThreshMapPath,8,1,0.3,{'thresh',[4,-3]}); %as above, but threshold INITIAL map using a positive threshold of >4 and negative <-3.
%        [ReThreshMapPath,ClusterMapPath] = ClusterReThresh(ThreshMapPath,8,1,0.3,{'thresh',[4,-3]},{'CLsize',18}); %as above, ALSO restruct clusters to >18 voxels otherwise delete + threshold INITIAL map using a positive threshold of >4 and negative <-3.
%        [ReThreshMapPath,ClusterMapPath] = ClusterReThresh(ThreshMapPath,8,1,0.3,{'CLsize',18}); %as above, but ONLY restrict clusters to >18 voxels otherwise delete.
%
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (20.02.2016): initial implementation based on previous tests

%% Check inputs
try
    if(ischar(ThreshMapPath))
        ThreshMapPath = cellstr(ThreshMapPath);
    else
        if(~iscellstr(ThreshMapPath))
            error('ThreshMapPath has to be string or cellstring!');
        end
    end
    TestPath = ThreshMapPath{1};
    [TestBaseDir,TestFName,TestExt] = fileparts(TestPath);
    TestExt = regexprep(TestExt,',\d*','');
    TestPath= [TestBaseDir,filesep,TestFName,TestExt];
    if(~exist(TestPath,'file'))
        error(['Can not find "',TestPath,'"!']);
    end
    clear TestPath TestBaseDir TestFName TestExt
catch CATCH_ThreshMapPath
    disp_catch(CATCH_ThreshMapPath,'ClusterReThresh>InputFile','CATCH_ThreshMapPath');
    ThreshMapPath = cellstr(spm_select(1,'image','Select thresholded map that should be rethresholded...'));
end

if(SearchDist_mm<=0)
    error('SearchDist_mm must be greater zero! (Distance in mm)');
end

if(UseConnectMat~=0)
    UseConnectMat = 1;
    disp('Restricting clustering of local maxima to connected voxels.');
end

if(RetainPercentage<=0)
    error('RetainPercentage must be greater zero (& smaller than 1)! (Percentage of top voxels to keep.)');
else
    if(RetainPercentage>=1)
        disp('WARNING: All voxels are retained after clustering!?');
    end        
end

switch(RetainMode)
    case {'Thresh','T','t','Threshold','threshold','thresh'}
        RetainMode = 'Thresh';
    case {'Voxels','voxels','Vox','vox','V','v'}
        RetainMode = 'Voxels';
    otherwise
        error(['RetainMode= "',RetainMode,'" is unknown.']);
end

%threshold initial map?
if(nargin>4)
    [ThreshInit, ClusterSizeThresh] = CheckXtaInputs(varargin);
else
    disp('Assuming input map is thresholded and that the cluster size limit it zero, i.e. taking all voxels that are non-zero.');
    ClusterSizeThresh = 0;
    ThreshInit        = [0 0];
end

%% load inputs
disp('Loading input map');
[XYZmm,StatsVals,XYZvox,Vol] = GetDataFromMap(ThreshMapPath{1},ThreshInit);

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
        ConnectMat{1}   = FindConnectedVoxels(XYZvox{1});
        disp('Reclustering connected voxels...');
        LocMaxStruct{1} = FindAllLocMax(XYZmm{1},StatsVals{1},SearchDist_mm,{'VoxelConnections',ConnectMat{1}});
        disp(' ');
    else
        ConnectMat{1}   = [];
        LocMaxStruct{1} = [];
    end
    %negative stats vals
    if(~isempty(XYZvox{2}))
        disp('Negative clusters will be examined');
        disp('Determining connectedness...');
        ConnectMat{2}   = FindConnectedVoxels(XYZvox{2});
        disp('Reclustering connected voxels...');
        LocMaxStruct{2} = FindAllLocMax(XYZmm{2},StatsVals{2},SearchDist_mm,{'VoxelConnections',ConnectMat{2}});
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
        disp(' ');
    else
        LocMaxStruct{1} = [];
    end
    %negative stats vals
    if(~isempty(XYZvox{2}))
        disp('Negative clusters will be examined');
        disp('Clustering voxels...');
        LocMaxStruct{2} = FindAllLocMax(XYZmm{2},StatsVals{2},SearchDist_mm);
        disp(' ');
    else
        LocMaxStruct{2} = [];
    end
end

%% apply ClusterSizeThresh?
if(ClusterSizeThresh>=1) %it is enough to have it >0 because a cluster is at least of size 1.
    XYZmmORG     = XYZmm;    %save a copy for comparison
    StatsValsORG = StatsVals;%save a copy for comparison
    XYZvoxORG    = XYZvox;   %save a copy for comparison
    for Ind = 1:2
        if(isstruct(LocMaxStruct{Ind})&&~isempty(LocMaxStruct{Ind}))
            LocMaxClusterSize = LocMaxStruct{Ind}.LocMaxClusterSize;
            if(any(LocMaxClusterSize<=ClusterSizeThresh))
                ClusterIndsDelete = sort(find(LocMaxClusterSize<=ClusterSizeThresh),'descend'); %when deleting entries we have to start with the last ones to preserve the structure of the previous ones, otherwise deletion will go wrong.
                if(Ind==1)
                    disp(['Throwing out ',num2str(length(find(LocMaxClusterSize<=ClusterSizeThresh))),' POSITIVE clusters (',num2str(ClusterIndsDelete(:)'),') of size <=',num2str(ClusterSizeThresh),'.']);
                else
                    disp(['Throwing out ',num2str(length(find(LocMaxClusterSize<=ClusterSizeThresh))),' NEGATIVE clusters (',num2str(ClusterIndsDelete(:)'),') of size <=',num2str(ClusterSizeThresh),'.']);
                end%if distinguishing positive & negative stats
                
                ClusterNo         = LocMaxStruct{Ind}.ClusterNo;        %(NVoxelx1)      The (final) cluster number per input voxel.
                LocMaxCoords_mm   = LocMaxStruct{Ind}.LocMaxCoords_mm;  %(NClustersx3)   The coordinate values of the local maxima.
                LocMaxClusterSize = LocMaxStruct{Ind}.LocMaxClusterSize;%(NClustersx1)   The number of voxels belonging to each cluster associated with each local maxima.
                CoG_Coords_mm     = LocMaxStruct{Ind}.CoG_Coords_mm;    %(NClustersx3)   The coordinate values of the CENTER of GRAVITY of the cluster, i.e. the average of all coordinates for the cluster of the local maxima.
                LocMaxStats       = LocMaxStruct{Ind}.LocMaxStats;      %(NClustersx3)   The first column is the statistics value at the local maxima and the second column is the average statistics values of the cluster and the third is the stdev.
                for IndCL = 1:length(ClusterIndsDelete)
                    ClusterNo(ClusterNo==ClusterIndsDelete(IndCL))= 0; %mark for deletion
                    LocMaxCoords_mm(ClusterIndsDelete(IndCL),:)   = []; %delete
                    LocMaxClusterSize(ClusterIndsDelete(IndCL))   = []; %delete
                    CoG_Coords_mm(ClusterIndsDelete(IndCL),:)     = []; %delete
                    LocMaxStats(ClusterIndsDelete(IndCL),:)       = []; %delete
                end
                %delete those that need to be thrown out
                XYZmm{Ind}(ClusterNo==0,:)   = []; %don't need these coords (mm)  any more
                StatsVals{Ind}(ClusterNo==0) = []; %don't need these stats values any more
                XYZvox{Ind}(ClusterNo==0,:)  = []; %don't need these voxel-coords any more
                ClusterNo(ClusterNo==0) = []; %delete
                %reassign clusters without gaps
                remainingClusters = unique(ClusterNo);
                for NewClusterInd = 1:length(unique(ClusterNo)) 
                    ClusterNo(ClusterNo==remainingClusters(NewClusterInd)) = NewClusterInd;
                    %disp(['OrgCL=',num2str(remainingClusters(NewClusterInd)),'-->',num2str(NewClusterInd),'=NewCL']); %for debug
                end
                
                %reinsert
                LocMaxStruct{Ind}.ClusterNo        = ClusterNo;        %(NVoxelx1)      The (final) cluster number per input voxel.
                LocMaxStruct{Ind}.LocMaxCoords_mm  = LocMaxCoords_mm;  %(NClustersx3)   The coordinate values of the local maxima.
                LocMaxStruct{Ind}.LocMaxClusterSize= LocMaxClusterSize;%(NClustersx1)   The number of voxels belonging to each cluster associated with each local maxima.
                LocMaxStruct{Ind}.CoG_Coords_mm    = CoG_Coords_mm;    %(NClustersx3)   The coordinate values of the CENTER of GRAVITY of the cluster, i.e. the average of all coordinates for the cluster of the local maxima.
                LocMaxStruct{Ind}.LocMaxStats      = LocMaxStats;      %(NClustersx3)   The first column is the statistics value at the local maxima and the second column is the average statistics values of the cluster.
            end%if throwing cluster that are too small
        end%if checking LocMaxStruct
    end%for-loop
end                

%% do rethresholding of each cluster: set all stats vals below top "RetainPercentage"% to zero --> then throw all such voxels out in the end and repeat stats. (NB: use ceil to determine number of voxels to keep.)
LocMaxStructFinal = LocMaxStruct;
XYZmmFinal     = XYZmm;    %save a copy for comparison AND for the output of the clusters BEFORE RETHRESHOLDING
StatsValsFinal = StatsVals;%save a copy for comparison AND for the output of the clusters BEFORE RETHRESHOLDING
XYZvoxFinal    = XYZvox;   %save a copy for comparison AND for the output of the clusters BEFORE RETHRESHOLDING
for Ind = 1:2
    if(isstruct(LocMaxStruct{Ind})&&~isempty(LocMaxStruct{Ind}))
        if(Ind==1) %positive stats vals
            disp(['Rethresholding POSITIVE clusters by retaining TOP ',num2str(RetainPercentage*100),'% relative to maximum.']);
            StatsSign =  1; %use this later to make calculations simpler, i.e. only coded once.
        else %negative stats vals
            disp(['Rethresholding NEGATIVE clusters by retaining "TOP" ',num2str(RetainPercentage*100),'% relative to maximum']);
            StatsSign = -1; %use this later to make calculations simpler, i.e. only coded once.
        end%if checking LocMaxStruct.
        
        ClusterNo = LocMaxStructFinal{Ind}.ClusterNo; %(NVoxelx1)      The (final) cluster number per input voxel.
        for IndCL = 1:max(ClusterNo)
            CurrCLstats    = StatsSign.*StatsValsFinal{Ind}(ClusterNo==IndCL); %should be positive due to StatsSign.
            IndsOfInterest = find(ClusterNo==IndCL); %the relevant indices %should be positive due to StatsSign.
            [SortCurrCLstats,SortInds] = sort(CurrCLstats,'descend');
            switch(RetainMode)
                case {'Thresh','T','t','Threshold','threshold','thresh'}
                    CurrThresh = SortCurrCLstats(max([1; floor(length(SortCurrCLstats)*(1-RetainPercentage))])); %using max([1; Ind]) ensures that we take a valid index at any configuration, i.e. at least the first/least significant.
                    
                    %replace voxels below current threshold with zero.
                    StatsValsFinal{Ind}((ClusterNo==IndCL)&((StatsSign.*StatsValsFinal{Ind})<CurrThresh)) = 0; %StatsSign makes this correct for pos/neg stats vals (logic: for those voxels to zero that are of the current cluster and "below" the current threshold)
                case {'Voxels','voxels','Vox','vox','V','v'}
                    ThrowInds = IndsOfInterest(SortInds(ceil(RetainPercentage*length(SortInds)):end));
                    StatsValsFinal{Ind}(ThrowInds) = 0; %set these to zero i.e. mark for deletion
                otherwise
                    error(['RetainMode= "',RetainMode,'" is unknown.']);
            end
        end%for-loop CLUSTERS
        %remove those that are zero.
        XYZmmFinal{Ind}( StatsValsFinal{Ind}==0,:)  = [];
        XYZvoxFinal{Ind}(StatsValsFinal{Ind}==0,:)  = [];
        LocMaxStructFinal{Ind}.ClusterNo(StatsValsFinal{Ind}==0) = [];
        StatsValsFinal{Ind}(StatsValsFinal{Ind}==0) = [];
        
        %redo calculations
        LocMaxStructFinal{Ind}.LocMaxClusterSize= zeros(max(LocMaxStructFinal{Ind}.ClusterNo),1); %(NClustersx1)   The number of voxels belonging to each cluster associated with each local maxima.
        LocMaxStructFinal{Ind}.CoG_Coords_mm    = zeros(max(LocMaxStructFinal{Ind}.ClusterNo),3); %(NClustersx3)   The coordinate values of the CENTER of GRAVITY of the cluster, i.e. the average of all coordinates for the cluster of the local maxima.
        LocMaxStructFinal{Ind}.LocMaxStats      = zeros(max(LocMaxStructFinal{Ind}.ClusterNo),3); %(NClustersx3)   The first column is the statistics value at the local maxima and the second column is the average statistics values of the cluster and the third is the stdev.
        for IndCL = 1:max(LocMaxStructFinal{Ind}.ClusterNo)
            LocMaxStructFinal{Ind}.LocMaxClusterSize(IndCL) = length(find(LocMaxStructFinal{Ind}.ClusterNo==IndCL));
            LocMaxStructFinal{Ind}.CoG_Coords_mm(IndCL,:)   = mean(XYZmmFinal{Ind}(LocMaxStructFinal{Ind}.ClusterNo==IndCL,:));
            LocMaxStructFinal{Ind}.LocMaxStats(IndCL,1)     = max( StatsValsFinal{Ind}(LocMaxStructFinal{Ind}.ClusterNo==IndCL));
            LocMaxStructFinal{Ind}.LocMaxStats(IndCL,2)     = mean(StatsValsFinal{Ind}(LocMaxStructFinal{Ind}.ClusterNo==IndCL));
            LocMaxStructFinal{Ind}.LocMaxStats(IndCL,3)     = std( StatsValsFinal{Ind}(LocMaxStructFinal{Ind}.ClusterNo==IndCL));
        end%for-loop CLUSTERS reassign
    end%if checking LocMaxStruct
end%for-loop Pos/Neg Stats

%% write out results to NIFTI-files
oVol = Vol;

[BaseDir,FNameOrg,ExtOrg] = fileparts(ThreshMapPath{1});
if(UseConnectMat)
    ConnStr = 'Conn';
else
    ConnStr = '';
end
AllVoxCoords  = []; %collect coords for output of pos & neg stats AFTER rethresholding
AllVoxStats   = []; %collect stats  for output of pos & neg stats AFTER rethresholding
ClusterMapPath= cell(2,1);
for Ind = 1:2
    if(isstruct(LocMaxStruct{Ind})&&~isempty(LocMaxStruct{Ind}))
        if(Ind==1)
            disp('Writing out POSITIVE STATS CLUSTERS BEFORE rethresholding.');
            StatsStr = 'PosStats';
        else
            disp('Writing out NEGATIVE STATS CLUSTERS BEFORE rethresholding.');
            StatsStr = 'NegStats';
        end
        if(ClusterSizeThresh>=1)
            CLsizeStr = ['CLsize',num2str(ClusterSizeThresh)];
        else
            CLsizeStr = '';
        end
        %Clusters are separated by stats sign (Pos/Neg)
        OutputPath          = [BaseDir,filesep,'C',num2str(SearchDist_mm),'mm',ConnStr,CLsizeStr,StatsStr,'_',FNameOrg,'.nii'];
        ClusterMapPath{Ind} = LocMaxClusters2NIFTI(LocMaxStruct{Ind}.ClusterNo,XYZvox{Ind},oVol,OutputPath);
        
        %collect voxels and stats for output (this one should go together, therefore collect).
        AllVoxCoords = [AllVoxCoords; XYZvoxFinal{Ind}];
        AllVoxStats  = [AllVoxStats;  StatsValsFinal{Ind}];
    end
end
OutputPath      = [BaseDir,filesep,'ReThresh',num2str(RetainPercentage*100),RetainMode,'C',num2str(SearchDist_mm),'mm',ConnStr,CLsizeStr,'_',FNameOrg,'.nii'];
ReThreshMapPath = LocMaxClusters2NIFTI(AllVoxStats,AllVoxCoords,oVol,OutputPath); %don't worry, it says LocMaxCluster2NIFTI, but it can deal with any kind of numbers just as long as we have the VOXEL-coords & a output volume struct (& a path&name).

%% done
disp(['Done rethresholding using clusters (search dist = ',num2str(SearchDist_mm),'mm',ConnMatStr,' and rethresholding Top ',num2str(RetainPercentage*100),'% of voxels per cluster).']);
disp('DONE.');
disp(' ');

end

%% subfunctions
%% CheckXtaInputs
function [ThreshInit, ClusterSizeThresh] = CheckXtaInputs(XtraInputs)
% This function check the extra inputs
ThreshInit       = [0 0];
ClusterSizeThresh= 0;

if(iscell(XtraInputs)) %should be becaue varargin should be a cell.
    for Ind = 1:length(XtraInputs)
        XtraInput = XtraInputs{Ind};
        if(iscell(XtraInput))
            if(length(XtraInput)==2)
                switch(XtraInput{1})
                    case 'thresh'
                        ThreshInit = XtraInput{2};
                        disp(['Will threshold INITIAL MAP with [',num2str(ThreshInit(1)),',',num2str(ThreshInit(2)),'].']);
                    case 'CLsize'
                        ClusterSizeThresh = XtraInput{2};
                        disp(['Will threshold CLUSTERS, -throwing out clusters with size <=',num2str(ClusterSizeThresh),'.']);
                    otherwise
                        error(['Command string unknown "',XtraInput{1},'"???']);
                end
            else
                error('Extra inputs must be cell-vectors of length 2!');
            end
        end
    end
else
    error('Extra inputs are not a cell-vector!!!???');
end

end

%% GetDataFromMap
function [XYZmm,StatsVals,XYZvox,Vol] = GetDataFromMap(ThreshMapPath,Thresh)
% This function loads the data from the NIFTI file, positive voxels are in the first cell and
% negative in the second entry of the cell.

%% get volume struct
Vol = spm_vol(ThreshMapPath);

%% extract data 
Data3D = Vol.private.dat(:,:,:,Vol.n(1));

%% assign data accordingly
XYZmm     = cell(2,1);
StatsVals = cell(2,1);
XYZvox    = cell(2,1);

%positive stats vals
if(any(Data3D(:)>Thresh(1)))
    StatsVals{1} = Data3D(Data3D(:)>Thresh(1)); %statistics values POSITIVE
    IndsStatsVals= find(Data3D(:)>Thresh(1));   %corresponding indices
    XYZvox{1}    = zeros(length(IndsStatsVals),3); %nVoxels-x-3 (spatial dimensions)
    [XYZvox{1}(:,1),XYZvox{1}(:,2),XYZvox{1}(:,3)] = ind2sub(Vol.dim,IndsStatsVals); %indices to 3D/subscript
    %do trafo from voxels to world (mm-coords)
    TestDiff = Vol.mat(1:3,1:3)-diag(diag(Vol.mat(1:3,1:3))); %if this contains only zeros then we can do this very fast!
    if(all(TestDiff(:)==0)) %rotation&scaling submatrix is diagonal --> we can do it fast with bsxfun
        XYZmm{1}=bsxfun(@plus,bsxfun(@times,XYZvox{1},diag(Vol.mat(1:3,1:3))'),Vol.mat(1:3,4)');
    else %not diagonal
        disp('WARNING(pos): VOXEL2WORLD 3D submatrix is NOT diagonal!!! (NB: usually any normalized (MNI) space should be diagonal in the 3D subspace.)'); %usually any normalized (MNI) space should be diagonal in the 3D subspace.
        XYZmm{1}=zeros(length(IndsStatsVals),3);
        for k = 1:length(IndsStatsVals)
            XYZmm{1}(k,:) = XYZvox{1}(k,:)*Vol.mat(1:3,1:3)+Vol.mat(1:3,4)';
        end
    end
end

%negative stats vals
if(any(Data3D(:)<Thresh(2)))
    StatsVals{2} = Data3D(Data3D(:)<Thresh(2)); %statistics values NEGATIVE
    IndsStatsVals= find(Data3D(:)<Thresh(2));   %corresponding indices
    XYZvox{2}    = zeros(length(IndsStatsVals),3); %nVoxels-x-3 (spatial dimensions)
    [XYZvox{2}(:,1),XYZvox{2}(:,2),XYZvox{2}(:,3)] = ind2sub(Vol.dim,IndsStatsVals); %indices to 3D/subscript
    %do trafo from voxels to world (mm-coords)
    TestDiff = Vol.mat(1:3,1:3)-diag(diag(Vol.mat(1:3,1:3))); %if this contains only zeros then we can do this very fast!
    if(all(TestDiff(:)==0)) %rotation&scaling submatrix is diagonal --> we can do it fast with bsxfun
        XYZmm{2}=bsxfun(@plus,bsxfun(@times,XYZvox{2},diag(Vol.mat(1:3,1:3))'),Vol.mat(1:3,4)'); %fast affine transformation using bsxfun (3D submatrix is diagonal)
    else %not diagonal
        disp('WARNING(neg): VOXEL2WORLD 3D submatrix is NOT diagonal!!! (NB: usually any normalized (MNI) space should be diagonal in the 3D subspace.)'); %usually any normalized (MNI) space should be diagonal in the 3D subspace.
        XYZmm{2}=zeros(length(IndsStatsVals),3);
        for k = 1:length(IndsStatsVals)
            XYZmm{2}(k,:) = XYZvox{2}(k,:)*Vol.mat(1:3,1:3)+Vol.mat(1:3,4)';
        end
    end
end

end