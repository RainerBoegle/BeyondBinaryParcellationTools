function [NormClMapPaths,ClusterMapPaths,VolNormCl] = CreateNormClustersMapFromStats(ThreshMapsList,ClusterSearchDist_mm,UseConnectMat,MinClusterSize)
% This function can be used to create NORMALIZED-Clusters Maps from Statistics Maps.
% The statistics maps are analysed for local maxima and associated clusters.
% Each voxels for each cluster is "weighted" based on the relative "statistical significance" to the local maxima.
% This NORMALIZES each cluster to have a value of "1" at the local maximum and adjusts the rest of each cluster relative to the local maxima.
% I.e. clusters with extremely high local maxima usually "bleed out" over neighboring regions,
% -in other words they loose their localization, and the normalization accounts for some of these
% spreads. For compact clusters the shape should still be maintained because they do not deviate
% much from the local maximum value.
% The end effect is that spread out clusters get more compact (as they should be) and compact
% clusters do not get penalized very much at the same time while equalizing the local maxima, 
% as all Clusters are significant and absolute value thresholding is not appropriate, see
% simulations (to do...).
% 
%
% For each StatsMap find initial clusters based on connected voxels. This forming of initial clusters
% is then proceeded by finding the local maxima withing these inital clusters and dividing the initial clusters
% into smaller clusters associated with the local maxima.
% NB: each cluster is tested for its size using "MinClusterSize" and clusters smaller than "MinClusterSize" will be thrown out.
% 
% For each such cluster find the statistics value of the local maximum and
% divide all StatsValues OF THE CURRENT CLUSTER by the value of the local maximum. 
% NB: for the StatsValues we first perform SignStats = sign(StatsValues); 
%     and use this to separate positive & negative clusters.
% % if(all(StatsVals>0))
% %     StatsSign = 1; %positive cluster
% % else
% %     StatsSign = -1; %negative cluster
% % end
% % StatsVals = StatsSign.*StatsVals;
% %
% %Result:
% % NormCl = StatsSign.*(StatsVals./max(StatsVals));
% 
%
%
%Inputs:
%    ThreshMapsList   cellstr (NImg-x-1) Paths to input stats images. (Should be thresholded already!)
%    ClSearchDist_mm  double   (1-x-1)   Distance for the creation of clusters, i.e. the distance that each voxel is allowed to search
%    (DEFAULT==8[mm])                    for another voxel that has a higher statistics value, IN EACH ITERATION. 
%         or                             If all voxels can not find any higher voxels in that distance  
%    (DEFAULT==2*Res[mm])                any more for several iterations (def==3), then convergence is reached,
%                                        i.e. all clusers have been formed. 
%                                        This procedure is done after all voxels have been pre-clustered
%                                        into cluters of connected voxels, that are then broken down into
%                                        smaller clusters (is possible) with the above approach. 
%                                        See help of functions ClusterThreshMap.m & FindAllLocMax.m.
%    UseConnectMat    boolean  (1-x-1)   Use connected voxels to form initial clusters. 
%                                        It is suggested to do this, except if the user wants to allow cluster formation
%                                        even from disconnected voxels.(NOT SUGGESTED.)
%    MinClusterSize   integer  (1-x-1)   The minimum cluster size for excluding clusters.
%
%
%Usage:
%      [NormClMapPaths,ClusterMapPaths,VolNormCl] = CreateNormClustersMapFromStats(ThreshMapsList,ClSearchDist_mm,UseConnectMat,MinClusterSize);
%      [NormClMapPaths,ClusterMapPaths,VolNormCl] = CreateNormClustersMapFromStats(ThreshMapsList,8,1,6);
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (16.03.2016): initial implementation

%% check inputs
try
    if(~iscellstr(ThreshMapsList))
        if(ischar(ThreshMapsList))
            ThreshMapsList = cellstr(ThreshMapsList);
        else
            error('"ThreshMapsList" must either be a char/string-array or cellstring!');
        end
    end
catch CATCH_ThreshMapsList
    disp_catch(CATCH_ThreshMapsList,'CreateNormClustersMapFromStats>Input:ThreshMapsList','CATCH_ThreshMapsList');
    ThreshMapsList = cellstr(spm_select(Inf,'image','Select thresholded maps for creation of RankMaps...'));
end

try
    if(ClusterSearchDist_mm<=0)
        disp('"ClusterSearchDist" must be greater than zero!');
    end
catch CATCH_ClusterSearchDist_mm
    disp_catch(CATCH_ClusterSearchDist_mm,'CreateNormClustersMapFromStats>Input:ClusterSearchDist_mm','CATCH_ClusterSearchDist_mm');
    Choice_ClusterSearchDist_mm = questdlg('Input "ClusterSearchDist_mm" manually or use 2*max(voxel-size)?[of first input image]','ClusterSearchDist_mm?','Input','2*max(Voxel-Size)','Input');
    switch(Choice_ClusterSearchDist_mm)
        case 'Input'
            answer_ClusterSearchDist_mm = inputdlg({'ClusterSearchDist[mm]= '},'ClusterSearchDist_mm?',1,{'8'});
            ClusterSearchDist_mm = str2num(answer_ClusterSearchDist_mm{1});
        otherwise
            NIItmp = nifti(ThreshMapsList{1});
            ClusterSearchDist_mm = 2*max(abs(diag(NIItmp.mat(1:3,1:3))));
    end
end

try
    if(~exist('UseConnectMat','var'))
        UseConnectMat = 1;
    else
        if(isempty(UseConnectMat))
            error('UseConnectMat has to be "0" or "1"!');
        else
            if(UseConnectMat<0)
                error('"UseConnectMat" must be "0" or "1"!');
            else
                if(UseConnectMat>0&&UseConnectMat<1)
                    UseConnectMat = 1;
                    disp('Setting UseConnectMat = 1;');
                elseif(UseConnectMat>1)
                    disp('Will set UseConnectMat = 1;');
                    UseConnectMat = 1;
                end
            end
        end
    end
catch CATCH_UseConnectMat
    disp_catch(CATCH_UseConnectMat,'CreateNormClustersMapFromStats>Input:UseConnectMat','CATCH_UseConnectMat');
    disp('Using default for UseConnectMat = 1; i.e. clustering voxels first by connectedness.');
    UseConnectMat = 1;
end

try
    if(~exist('MinClusterSize','var'))
        MinClusterSize = 0; %do not use.
    else
        if(isempty(MinClusterSize))
            MinClusterSize = 0; %do not use.
        else
            if(MinClusterSize<0)
                error('"MinClusterSize" can not be smaller than zero!');
            else
                disp(['Using "MinClusterSize"= ',num2str(MinClusterSize),'.']);
            end
        end
    end
catch CATCH_MinClusterSize
    disp_catch(CATCH_MinClusterSize,'CreateNormClustersMapFromStats>Input:MinClusterSize','CATCH_MinClusterSize');
    MinClusterSize = 0;
end
if(rem(MinClusterSize,1)~=0)
    disp(['"MinClusterSize" (',num2str(MinClusterSize),') must be an integer!']);
    MinClusterSize = round(MinClusterSize);
    disp(['Will round it, i.e. MinClusterSize= ',num2str(MinClusterSize),'.']);
end

%% cluster input ThreshMaps AND create RankMaps
Thresh = [0 0]; %always assume that maps are thresholded already

NormClMapPaths  = cell(length(ThreshMapsList),1);
ClusterMapPaths = cell(length(ThreshMapsList),1);
VolNormCl       = cell(length(ThreshMapsList),1);
disp(' ');
for IndMap = 1:length(ThreshMapsList)
    [BaseDir,FName,Ext] = fileparts(ThreshMapsList{IndMap});
    disp(' ');
    disp(['Creating NormClMap(',num2str(IndMap,['%0',num2str(max([ceil(log10(length(ThreshMapsList))); 2])),'d']),'of',num2str(length(ThreshMapsList),['%0',num2str(max([ceil(log10(length(ThreshMapsList))); 2])),'d']),') for Input ',num2str(IndMap),': "',FName,Ext,'" (in directory "',BaseDir,'")...']);
    %Cluster
    [ClusterMapPaths{IndMap},VolCL,LocMaxStruct,XYZmm,StatsVals,XYZvox] = ClusterThreshMap(ThreshMapsList{IndMap},Thresh,UseConnectMat,ClusterSearchDist_mm,MinClusterSize);
    
    %Create NormCl Values for NormClMap
    [NormClTotal,NormClVoxTotal] = NormalizeAllClusters(LocMaxStruct,XYZmm,StatsVals,XYZvox);
    
    %write out NormClMap
    [BaseDir,FNameCL,Ext] = fileparts(VolCL{1}.fname);
    FNameCL = regexprep(regexprep(FNameCL,'PosStats',''),'NegStats','');
    disp(['Writing out Normalized Clusters "NormCl',FNameCL,Ext,'" to directory "',BaseDir,'"...']); 
    [NormClMapPaths{IndMap},VolNormCl{IndMap}] = LocMaxClusters2NIFTI(NormClTotal,NormClVoxTotal,VolCL{1},[BaseDir,filesep,'NormCl',FNameCL,Ext]);
end

%% Done.
disp('DONE.');
disp(' ');

end

%% subfunctions
%% MakeNormClForAllClusters
function [NormClTotal,NormClVoxTotal,NormClVoxTotal_mm] = NormalizeAllClusters(LocMaxStruct,XYZmm,StatsVals,XYZvox)
%This function will go over all clusters (positive and negative) and create the NormCl data and
%return it together with the appropriate Voxel- & mm-Coordinates to write these out to NIFTI.


%% calculations
NormCl       = cell(2,1); %positive & negative --> later combine them if not empty still.
NormClVox    = cell(2,1); %positive & negative --> later combine them if not empty still.
NormClVox_mm = cell(2,1); %positive & negative --> later combine them if not empty still.
for Ind = 1:2
    if(isstruct(LocMaxStruct{Ind})&&~isempty(LocMaxStruct{Ind}))
        ClusterNo = LocMaxStruct{Ind}.ClusterNo;
        CLInds    = unique(ClusterNo);
        if(Ind==1)
            disp('Positive clusters will be examined');
            for CLIndex = 1:length(CLInds)
                disp([num2str(CLIndex,['%0',num2str(max([ceil(log10(length(CLInds))); 2])),'d']),'of',num2str(length(CLInds),['%0',num2str(max([ceil(log10(length(CLInds))); 2])),'d'])]);
                CurrCL = CLInds(CLIndex);
                CurrCL_StatsVals = StatsVals{Ind}(ClusterNo==CurrCL);
                CurrCL_Vox_mm    =     XYZmm{Ind}(ClusterNo==CurrCL,:);
                CurrCL_Vox       =    XYZvox{Ind}(ClusterNo==CurrCL,:);
                NormCl{      Ind} = [   NormCl{   Ind}; CreateNormCl(CurrCL_StatsVals)];
                NormClVox{   Ind} = [NormClVox{   Ind}; CurrCL_Vox];
                NormClVox_mm{Ind} = [NormClVox_mm{Ind}; CurrCL_Vox_mm];
            end
        else
            disp('Negative clusters will be examined');
            for CLIndex = 1:length(CLInds)
                disp([num2str(CLIndex,['%0',num2str(max([ceil(log10(length(CLInds))); 2])),'d']),'of',num2str(length(CLInds),['%0',num2str(max([ceil(log10(length(CLInds))); 2])),'d'])]);
                CurrCL = CLInds(CLIndex);
                CurrCL_StatsVals = StatsVals{Ind}(ClusterNo==CurrCL);
                CurrCL_Vox_mm    =     XYZmm{Ind}(ClusterNo==CurrCL,:);
                CurrCL_Vox       =    XYZvox{Ind}(ClusterNo==CurrCL,:);
                NormCl{      Ind} = [   NormCl{   Ind}; CreateNormCl(CurrCL_StatsVals)];
                NormClVox{   Ind} = [NormClVox{   Ind}; CurrCL_Vox];
                NormClVox_mm{Ind} = [NormClVox_mm{Ind}; CurrCL_Vox_mm];
            end
        end
    else
        if(Ind==1)
            disp('NO positive clusters!');
        else
            disp('NO negative clusters!');
        end
    end
end

%% collect results %NormClTotal,NormClVoxTotal,NormClVoxTotal_mm
disp('Collecting results...');
NormClTotal      = [NormCl{      1}; NormCl{      2}];%init
NormClVoxTotal   = [NormClVox{   1}; NormClVox{   2}];%init
NormClVoxTotal_mm= [NormClVox_mm{1}; NormClVox_mm{2}];%init


%% check output
%% final check of size: NormClTotal
if(size(NormClTotal,1)==1&&size(NormClTotal,2)~=1)
    NormClTotal = NormClTotal';
elseif(size(NormClTotal,1)==1&&size(NormClTotal,2)==1)
    disp('Only one NormClTotal here???');
elseif(size(NormClTotal,1)~=1&&size(NormClTotal,2)==1)
    %correct
else
    size(NormClTotal)
    error('NormClTotal has strange size!');
end

%% final check of size: NormClVoxTotal
if(size(NormClVoxTotal,2)~=3&&size(NormClVoxTotal,1)==3)
    NormClVoxTotal = NormClVoxTotal';
elseif(size(NormClVoxTotal,1)==1&&size(NormClVoxTotal,2)==3)
    disp('Only one NormClVoxTotal here???');
elseif((size(NormClVoxTotal,1)~=1||size(NormClVoxTotal,1)~=3)&&size(NormClVoxTotal,2)==3)
    %correct
else
    size(NormClVoxTotal)
    error('NormClVoxTotal has strange size!');
end


        
end

%% CreateNormCl
function NormCl = CreateNormCl(StatsVals)
% if(all(StatsVals>0))
%     StatsSign = 1; %positive cluster
% else
%     StatsSign = -1; %negative cluster
% end
% StatsVals = StatsSign.*StatsVals;
%
%Result:
% NormCl = StatsSign.*(StatsVals./max(StatsVals));
%
%Inputs:
%       StatsVals   (NVox-x-1)  Vector containing the stats values of the current cluster of size NVox 


%% check inputs & get StatsSign
if(all(StatsVals>0))
    StatsSign = 1; %positive cluster
else
    StatsSign = -1; %negative cluster
end
if(any((StatsSign.*StatsVals)<0)) %this is not allowed (i.e. there is a mix of positive and negative values in this cluster) --> error
    error('Current cluster has positive and negative statistic values!!! This is not allowed, clusters should be either positive or negative, -not a mix!');
else
    StatsVals = StatsSign.*StatsVals; %replace stats vals --> easy use for later (and we do not return them or change the original ones therefore not problem here, -as long as this is not chenged in the program later...)
end

%% do normalization
disp(['Scaling stats for Cluster by local maximum value max(StatsVals)= ',num2str(max(StatsVals)),'(a.u.) [range(StatsVals)= ',num2str(range(StatsVals)),'(a.u.) --> [',num2str(min(StatsVals)),'...',num2str(max(StatsVals)),']<==>[',num2str(min(StatsVals)./max(StatsVals)),'...',num2str(max(StatsVals)./max(StatsVals)),']].']);
NormCl = StatsSign.*(StatsVals./max(StatsVals)); %inverse distance INCLUDING THE SIGN OF THE STATS!!!.
   

%% final check of size
if(size(NormCl,1)==1&&size(NormCl,2)~=1)
    NormCl = NormCl';
elseif(size(NormCl,1)==1&&size(NormCl,2)==1)
    disp('Only one NormCl here???');
elseif(size(NormCl,1)~=1&&size(NormCl,2)==1)
    %correct
else
    size(NormCl)
    error('NormCl has strange size!');
end

end