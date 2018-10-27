function [InvDistMapPaths,ClusterMapPaths,VolInvDist] = CreateInvDistMapFromStats(ThreshMapsList,ClusterSearchDist_mm,ExpectedClusterDistFromLocMax)
% This function can be used to create INVERSE-Distance Maps from Statistics Maps.
% The statistics maps are analysed for local maxima and associated clusters.
% Each voxels for each cluster is "weighted" based on the relative "distance" from the local maxima.
% The "distance" is the combination of the exponential of the scaled distance from the local maxima
% statistics value and local maxima coordinate, 
% i.e. D = sqrt(D1.*D2); 
% with D1 = (exp(StatsDistSc)-1)./(exp(1)-1);
%  &&  D2 = (exp(VoxDistSc)  -1)./(exp(1)-1); 
%
% and 
% StatsDist =     abs( StatsVals -        LocMaxStats                       );
% VoxDist   =sqrt(sum((Vox_mm    - repmat(LocMax_mm,size(Vox_mm,1),1)).^2,2));
%  
% with
% ScaleStats = weights(1).*range(StatsVals); %scale for stats of a cluster should be related to maximum difference from local maxima
% ScaleDist  = max([ExpClDistLocMax; (weights(2).*max(VoxDist))]); %scale for distance of voxels to the local maxima of a cluster should be related to the maximum distance from the local maxima but limited on the lower end by the expected distance.
% with weights(1) = 1/2 & weights(2) = abs(1-range(StatsVals)/max(StatsVals));
% 
% StatsDistSc = StatsDist./ScaleStats;   %Distance relative to scale (statistics distance)
% VoxDistSc   =   VoxDist./ScaleDist; %Distance relative to scale (coordinate distance)
% From this distance the INVERSE DISTANCE is created on the basis of InvDist = 1./(1+D).
%
% From this distance the INVERSE DISTANCE is created on the basis of InvDist = 1./(1+D).
% I.e. each cluster starts at 1 for the peak and then descends to 1/2 at the combined ScaleLength of the stats and distances.
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
%    ExpClDistLocMax  double  (1-x-1)    Distance from LocalMaxima that is expected for a cluster.
%    (DEFAULT==2*Res[mm])                This distance will be used as the scale length (ScaleDist see above),
%            or                          if the determined scale length is smaller than this.
%    (DEFAULT==6[mm])                    I.e. small clusters are not excluded but large clusters "bleeding out" 
%                                        from huge significant local maxima are diminished.
%
%Usage:
%      [InvDistMapPaths,ClusterMapPaths,VolInvDist] = CreateInvDistMapFromStats(ThreshMapsList,ClSearchDist_mm,ExpClDistLocMax);
%      [InvDistMapPaths,ClusterMapPaths,VolInvDist] = CreateInvDistMapFromStats(ThreshMapsList,8,6);
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
    disp_catch(CATCH_ThreshMapsList,'CreateInvDistMapFromStats>Input:ThreshMapsList','CATCH_ThreshMapsList');
    ThreshMapsList = cellstr(spm_select(Inf,'image','Select thresholded maps for creation of RankMaps...'));
end

try
    if(ClusterSearchDist_mm<=0)
        disp('"ClusterSearchDist" must be greater than zero!');
    end
catch CATCH_ClusterSearchDist_mm
    disp_catch(CATCH_ClusterSearchDist_mm,'CreateInvDistMapFromStats>Input:ClusterSearchDist_mm','CATCH_ClusterSearchDist_mm');
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
    if(~exist('ExpectedClusterDistFromLocMax','var'))
        NIItmp = nifti(ThreshMapsList{1});
        ExpectedClusterDistFromLocMax = 3*max(abs(diag(NIItmp.mat(1:3,1:3))));
    else
        if(isempty(ExpectedClusterDistFromLocMax))
            NIItmp = nifti(ThreshMapsList{1});
            ExpectedClusterDistFromLocMax = 3*max(abs(diag(NIItmp.mat(1:3,1:3))));
        else
            if(ExpectedClusterDistFromLocMax<=0)
                error('"ExpectedClusterDistFromLocMax" must be greater than zero and should best be 2*ImagingResolution in mm.');
            end
        end
    end
catch CATCH_ExpectedClusterDistFromLocMax
    disp_catch(CATCH_ExpectedClusterDistFromLocMax,'CreateInvDistMapFromStats>Input:ExpectedClusterDistFromLocMax','CATCH_ExpectedClusterDistFromLocMax');
    NIItmp = nifti(ThreshMapsList{1});
	ExpectedClusterDistFromLocMax = 3*max(abs(diag(NIItmp.mat(1:3,1:3))));
end

%% cluster input ThreshMaps AND create RankMaps
UseConnectMat = 1;     %always make connected clusters first.
Thresh        = [0 0]; %always assume that maps are thresholded already

InvDistMapPaths = cell(length(ThreshMapsList),1);
ClusterMapPaths = cell(length(ThreshMapsList),1);
VolInvDist      = cell(length(ThreshMapsList),1);
for IndMap = 1:length(ThreshMapsList)
    [BaseDir,FName,Ext] = fileparts(ThreshMapsList{IndMap});
    disp(['Creating InvDistMap(',num2str(IndMap,['%0',num2str(max([ceil(log10(length(ThreshMapsList))); 2])),'d']),'of',num2str(length(ThreshMapsList),['%0',num2str(max([ceil(log10(length(ThreshMapsList))); 2])),'d']),') for Input ',num2str(IndMap),': "',FName,Ext,'" (in directory "',BaseDir,'")...']);
    %Cluster
    [ClusterMapPaths{IndMap},VolCL,LocMaxStruct,XYZmm,StatsVals,XYZvox] = ClusterThreshMap(ThreshMapsList{IndMap},Thresh,UseConnectMat,ClusterSearchDist_mm);
    
    %Create InvDist Values for InvDistMap
    %OLD: weights = [1 1; 1 1]; %weigh median & mad for stats & distance equally; e.g. [1 0; 1 0] use only medians.
    weights = {1/2, 'abs(1-range/max)'}; %Alternative "NormalizeClusters" = 'Simple 1/max'; %OLD-Final = {1/2, 'abs(1-range/max)'}; %OLD4 = {1/2, '1/max'}; %OLD3 = {1/2, 'min/max'}; %OLD1 = [1/2 1/2]; %OLD2 = [1/2 1/5]; %or {1/2, 'min/max'}; %using weights as a cell and including 'min/max' will change the weights based on the min(StatsVals)/max(StatsVals). %the first entry is the weight for the stats distance based on the range of stats values i.e. weights(1).*range(StatVals); %the second entry is the weight for the spatial distance when calculating it, BUT keep in mind that this is limited on the lower end by ExpectedClusterDistFromLocMax.
    [InvDistTotal,InvDistVoxTotal] = MakeInvDistForAllClusters(LocMaxStruct,XYZmm,StatsVals,XYZvox,weights,ExpectedClusterDistFromLocMax);
    
    %write out InvDistMap
    [BaseDir,FNameCL,Ext] = fileparts(VolCL{1}.fname);
    FNameCL = regexprep(regexprep(FNameCL,'PosStats',''),'NegStats','');
    disp(['Writing out InvDistMap "InvDistMap_ExpCSiz',num2str(ExpectedClusterDistFromLocMax),'mm',FNameCL,Ext,'" to directory "',BaseDir,'"...']); 
    [InvDistMapPaths{IndMap},VolInvDist{IndMap}] = LocMaxClusters2NIFTI(InvDistTotal,InvDistVoxTotal,VolCL{1},[BaseDir,filesep,'InvDistMap_ExpCSiz',num2str(ExpectedClusterDistFromLocMax),'mm',FNameCL,Ext]);
end

%% Done.
disp('DONE.');
disp(' ');

end

%% subfunctions
%% MakeInvDistForAllClusters
function [InvDistTotal,InvDistVoxTotal,InvDistVoxTotal_mm] = MakeInvDistForAllClusters(LocMaxStruct,XYZmm,StatsVals,XYZvox,weights,ExpectedClusterDistFromLocMax)
%This function will go over all clusters (positive and negative) and create the InvDist data and
%return it together with the appropriate Voxel- & mm-Coordinates to write these out to NIFTI.

%% info for weigths
if(ischar(weights))
    if(strcmp(weights,'Simple 1/max'))
        disp('Will use simple algorithm to adjust clusters by the maximum/Local Maximum.');
    end
else
    if(iscell(weights))
        weightsInput = weights;
        %first entry
        if(ischar(weightsInput{1}))
            switch(weightsInput{1})
                case 'min/max'
                    disp(['Method-weights(1): "',weightsInput{1},'".']);
                    %can't calculate it here just info, later use formula: weights(1) = min(StatsVals)/max(StatsVals);
                case '1/max'
                    disp(['Method-weights(1): "',weightsInput{1},'".']);
                    %can't calculate it here just info, later use formula: weights(1) = 1/max(StatsVals);
                case 'abs(1-range/max)'
                    disp(['Method-weights(1): "',weightsInput{1},'".']);
                otherwise
                    disp('Method-weights(1): "fixed".');
                    %can't calculate it here just info, later use formula: weights(1) = str2num(weightsInput{1});
            end
        elseif(isscalar(weightsInput{1}))
            disp('Method-weights(1) is "fixed".');
            %can't calculate it here just info, later use formula: weights(1) = weightsInput{1};
        else
            error('"weights{1} is weird!!! Check the inputs in debugging mode.');
        end
        %second entry
        if(ischar(weightsInput{2}))
            switch(weightsInput{2})
                case 'min/max'
                    disp(['Method-weights(2): "',weightsInput{2},'".']);
                    %can't calculate it here just info, later use formula: weights(2) = min(StatsVals)/max(StatsVals);
                case '1/max'
                    disp(['Method-weights(2): "',weightsInput{2},'".']);
                    %can't calculate it here just info, later use formula: weights(2) = 1/max(StatsVals);
                case 'abs(1-range/max)'
                    disp(['Method-weights(2): "',weightsInput{2},'".']);
                otherwise
                    disp('Method-weights(2): "fixed".');
                    %can't calculate it here just info, later use formula: weights(2) = str2num(weightsInput{2});
            end
        elseif(isscalar(weightsInput{2}))
            disp('Method-weights(2) is "fixed".');
            %can't calculate it here just info, later use formula: weights(2) = weightsInput{2};
        else
            error('"weights{2} is weird!!! Check the inputs in debugging mode.');
        end
    end
end

%% calculations
InvDist       = cell(2,1); %positive & negative --> later combine them if not empty still.
InvDistVox    = cell(2,1); %positive & negative --> later combine them if not empty still.
InvDistVox_mm = cell(2,1); %positive & negative --> later combine them if not empty still.
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
                InvDist{      Ind} = [   InvDist{   Ind}; CreateInvDist(CurrCL_StatsVals,CurrCL_Vox_mm,weights,ExpectedClusterDistFromLocMax)];
                InvDistVox{   Ind} = [InvDistVox{   Ind}; CurrCL_Vox];
                InvDistVox_mm{Ind} = [InvDistVox_mm{Ind}; CurrCL_Vox_mm];
            end
        else
            disp('Negative clusters will be examined');
            for CLIndex = 1:length(CLInds)
                disp([num2str(CLIndex,['%0',num2str(max([ceil(log10(length(CLInds))); 2])),'d']),'of',num2str(length(CLInds),['%0',num2str(max([ceil(log10(length(CLInds))); 2])),'d'])]);
                CurrCL = CLInds(CLIndex);
                CurrCL_StatsVals = StatsVals{Ind}(ClusterNo==CurrCL);
                CurrCL_Vox_mm    =     XYZmm{Ind}(ClusterNo==CurrCL,:);
                CurrCL_Vox       =    XYZvox{Ind}(ClusterNo==CurrCL,:);
                InvDist{      Ind} = [   InvDist{   Ind}; CreateInvDist(CurrCL_StatsVals,CurrCL_Vox_mm,weights,ExpectedClusterDistFromLocMax)];
                InvDistVox{   Ind} = [InvDistVox{   Ind}; CurrCL_Vox];
                InvDistVox_mm{Ind} = [InvDistVox_mm{Ind}; CurrCL_Vox_mm];
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

%% collect results
disp('Collecting results...');
InvDistTotal      = [InvDist{      1}; InvDist{      2}];%init
InvDistVoxTotal   = [InvDistVox{   1}; InvDistVox{   2}];%init
InvDistVoxTotal_mm= [InvDistVox_mm{1}; InvDistVox_mm{2}];%init


%% check output
%% final check of size: InvDistTotal
if(size(InvDistTotal,1)==1&&size(InvDistTotal,2)~=1)
    InvDistTotal = InvDistTotal';
elseif(size(InvDistTotal,1)==1&&size(InvDistTotal,2)==1)
    disp('Only one InvDistTotal here???');
elseif(size(InvDistTotal,1)~=1&&size(InvDistTotal,2)==1)
    %correct
else
    size(InvDistTotal)
    error('InvDistTotal has strange size!');
end

%% final check of size: InvDistVoxTotal
if(size(InvDistVoxTotal,2)~=3&&size(InvDistVoxTotal,1)==3)
    InvDistVoxTotal = InvDistVoxTotal';
elseif(size(InvDistVoxTotal,1)==1&&size(InvDistVoxTotal,2)==3)
    disp('Only one InvDistVoxTotal here???');
elseif((size(InvDistVoxTotal,1)~=1||size(InvDistVoxTotal,1)~=3)&&size(InvDistVoxTotal,2)==3)
    %correct
else
    size(InvDistVoxTotal)
    error('InvDistVoxTotal has strange size!');
end


        
end

%% CreateInvDist
function InvDist = CreateInvDist(StatsVals,Vox_mm,weights,ExpectedClusterDistFromLocMax)
% i.e. D = sqrt(D1.*D2); 
% with D1 = (exp(StatsDistSc)-1)./(exp(1)-1);
%  &&  D2 = (exp(VoxDistSc)  -1)./(exp(1)-1); 
%
% and 
% if(all(StatsVals>0))
%     StatsSign = 1; %positive cluster
% else
%     StatsSign = -1; %negative cluster
% end
% StatsVals = StatsSign.*StatsVals;
%
% %% get stats dist and VoxDist
% StatsDist =     abs( StatsVals -        LocMaxStats                       );
% VoxDist   =sqrt(sum((Vox_mm    - repmat(LocMax_mm,size(Vox_mm,1),1)).^2,2));
%  
% %% create the scale length for the stats & distances of voxels.
% ScaleStats = weights(1).*range(StatsVals); %scale for stats of a cluster should be related to maximum difference from local maxima
% ScaleDist  = max([ExpectedClusterDistFromLocMax; (weights(2).*max(VoxDist))]); %scale for distance of voxels to the local maxima of a cluster should be related to the maximum distance from the local maxima but limited on the lower end by the expected distance.
% 
% %% apply the scale length
% StatsDistSc = StatsDist./ScaleStats;   %Distance relative to scale (statistics distance)
% VoxDistSc   =   VoxDist./ScaleDist; %Distance relative to scale (coordinate distance)
%
% From this distance the INVERSE DISTANCE is created on the basis of InvDist = StatsSign.*(1./(1+D)).
%
%Inputs:
%       StatsVals   (NVox-x-1)  Vector containing the stats values of the current cluster of size NVox 
%       Vox_mm      (NVox-x-3)  The XYZ-coordinates in mm (assuming MNI space but any space that is metric and has an equivalent to mm distance is also fine, -I guess)
%       weights      (2-x-2)    How much weight is assigned to median & mad for stats & distance; [1 1; 1 1]; means equal weights and e.g. [1 0; 1 0] means use only medians.

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

if(size(Vox_mm,1)~=length(StatsVals))
    error(['"Vox_mm"(',num2str(size(Vox_mm,1)),') does not contain the same number of voxel coordinates as stats values in "StatsVals"(',num2str(length(StatsVals)),')!']);
end


%% check weights
if(ischar(weights))
    if(strcmp(weights,'Simple 1/max'))
        %% create "inverse distance"
        InvDist = StatsSign.*(StatsVals./max(StatsVals)); %inverse distance INCLUDING THE SIGN OF THE STATS!!!.
    else
        error(['weights selection "',weights,'" unknown!']);
    end
else
    if(iscell(weights))
        weightsInput = weights; clear weights
        %first entry
        if(ischar(weightsInput{1}))
            switch(weightsInput{1})
                case 'min/max'
                    %No status message here, has been done before: disp(['Method-weights(1): "',weightsInput{1},'".']);
                    weights(1) = min(StatsVals)/max(StatsVals);
                case '1/max'
                    %No status message here, has been done before: disp(['Method-weights(1): "',weightsInput{1},'".']);
                    weights(1) = 1/max(StatsVals);
                case 'abs(1-range/max)'
                    weights(1) = abs(1-range(StatsVals)/max(StatsVals));
                otherwise
                    %No status message here, has been done before: disp('Method-weights(1): "fixed".');
                    weights(1) = str2num(weightsInput{1});
            end
        elseif(isscalar(weightsInput{1}))
            %No status message here, has been done before: disp('Method-weights(1) is "fixed".');
            weights(1) = weightsInput{1};
        else
            error('"weights{1} is weird!!! Check the inputs in debugging mode.');
        end
        %second entry
        if(ischar(weightsInput{2}))
            switch(weightsInput{2})
                case 'min/max'
                    %No status message here, has been done before: disp(['Method-weights(2): "',weightsInput{2},'".']);
                    weights(2) = min(StatsVals)/max(StatsVals);
                case '1/max'
                    %No status message here, has been done before: disp(['Method-weights(2): "',weightsInput{2},'".']);
                    weights(2) = 1/max(StatsVals);
                case 'abs(1-range/max)'
                    weights(2) = abs(1-range(StatsVals)/max(StatsVals));
                otherwise
                    %No status message here, has been done before: disp('Method-weights(2): "fixed".');
                    weights(2) = str2num(weightsInput{2});
            end
        elseif(isscalar(weightsInput{2}))
            %No status message here, has been done before: disp('Method-weights(2) is "fixed".');
            weights(2) = weightsInput{2};
        else
            error('"weights{2} is weird!!! Check the inputs in debugging mode.');
        end
    end
    
    %% get LocMaxStats and LocMax_mm
    [LocMaxStats,IndLocMax] = max(StatsVals);
    LocMax_mm              = Vox_mm(IndLocMax,:);
    
    %% get stats dist and VoxDist
    StatsDist =     abs( StatsVals -        LocMaxStats                       );
    VoxDist   =sqrt(sum((Vox_mm    - repmat(LocMax_mm,size(Vox_mm,1),1)).^2,2));
    
    %% create the scale length for the stats & distances of voxels.
    % ScaleStats =                             1/2.*(weights(1,1).*median(StatsDist(StatsVals>median(StatsVals)))+ weights(1,2).*mad(StatsDist(StatsVals>median(StatsVals))));
    % ScaleDist  = (min(StatsVals)/max(StatsVals)).*(weights(2,1).*median(  VoxDist(StatsVals>median(StatsVals)))+ weights(2,2).*mad(  VoxDist(StatsVals>median(StatsVals))));
    % %NB: still need to work on this but the intuition is that a huge range should lead to a huge cluster that should not be some huge and different in values,
    % %    i.e. punish this by making the scale proportionetely smaller, i.e. emphasize the deviations for this cluster more.
    % 2/3*min maybe???
    
    ScaleStats = weights(1).*range(StatsVals); %scale for stats of a cluster should be related to maximum difference from local maxima
    ScaleDist  = max([ExpectedClusterDistFromLocMax; (weights(2).*max(VoxDist))]); %scale for distance of voxels to the local maxima of a cluster should be related to the maximum distance from the local maxima but limited on the lower end by the expected distance.
    
    disp(['ScaleStats= ',num2str(ScaleStats),'(a.u.)=weights(1).*range(StatsVals); [max(StatsVals)=',num2str(max(StatsVals)),'(a.u.) & range(StatsVals)=',num2str(range(StatsVals)),'(a.u.)]. (weights(1)= ',num2str(weights(1)),')']);
    if(ScaleDist==ExpectedClusterDistFromLocMax)
        disp(['ScaleDist will be "ExpectedClusterDistFromLocMax"= ',num2str(ScaleDist),'mm. (Alternative would have been ',num2str((weights(2).*max(VoxDist))),'mm (=weights(2).*max(VoxDist)). weights(2)= ',num2str(weights(2)),' NOT USED!)']);
    else
        disp(['ScaleDist= ',num2str((weights(2).*max(VoxDist))),'mm (=weights(2).*max(VoxDist)). (weights(2)= ',num2str(weights(2)),')']);
    end
    
    
    %% apply the scale length
    StatsDistSc = StatsDist./ScaleStats;   %Distance relative to scale (statistics distance)
    VoxDistSc   =   VoxDist./ScaleDist; %Distance relative to scale (coordinate distance)
    
    %% create Transformed distances
    D1 = (exp(StatsDistSc)-1)./(exp(1)-1); %transform distance exponential; start at zero and be at 1 at scale distance then just grow from there exponentially
    D2 = (exp(VoxDistSc)  -1)./(exp(1)-1); %transform distance exponential; start at zero and be at 1 at scale distance then just grow from there exponentially
    
    %% multiply TRANSFORMED distances
    D = sqrt(D1.*D2); %product distance as square root (geometric mean)
    
    %% create inverse distance
    InvDist = StatsSign.*(1./(1+D)); %inverse distance INCLUDING THE SIGN OF THE STATS!!!.
end

%% final check of size
if(size(InvDist,1)==1&&size(InvDist,2)~=1)
    InvDist = InvDist';
elseif(size(InvDist,1)==1&&size(InvDist,2)==1)
    disp('Only one InvDist here???');
elseif(size(InvDist,1)~=1&&size(InvDist,2)==1)
    %correct
else
    size(InvDist)
    error('InvDist has strange size!');
end

end