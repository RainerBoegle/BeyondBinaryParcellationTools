function [LocMaxStruct,MapExtractStruct,ClusterNIIOutputPath] = ClusterStatsMap(StatsMapPath,ThresholdSettings,EffectiveDist_mm,UseConnectMat,ClusterNIIOutputDir)
% This function provides the interface for all tools that are needed to determine all local maxima
% of a thresholded statistics map (something like a z-, t- or F-test values) and assigning all voxels
% around those local maxima to form clusters.
% The clustering can be done on the basis of a NIFTI-file or a SPM.mat & Contrast-Settings.
%
%INPUTS:
%  StatsMapPath       <--   Can either be the location (string) of a SPM.mat-file or a "StatsMap".nii/img&hdr-file.
%                             NB: THE TYPE OF INPUT DETERMINES WHAT IS EXPECTED IN "ThresholdSettings". 
%  ThresholdSettings  <--   In case of a NIFTI-file location in "StatsMapPath" this should be a (2-x-1) vector
%                           with the negative and positive thresholds. In case one threshold does not apply, replace it with "nan". 
%                           E.g. ThresholdSettings= [-4.496,4.496] OR ThresholdSettings= [nan,4.496].
%
%                           In case of a SPM.mat "ThresholdSettings" has to be struct, containing the following fields:
%                               .Contrast     <-- For selecting the contrast. This should either be a string or a number. 
%                                                 If string then it should fit the name of the desired contrast. 
%                                                 If number then it indicates the contrast via this number.'
%                                                     NB: If it is left empty (i.e. [] or missing), then you will be asked to select it manually
%                               .MultTestCor  <-- The multiple correction method to be used. Should be a string, either 'FWE' or 'none'. 
%                                                 %"FDR" might work, but I doubt it, -not checked! May need an edit of spm_defaults.m "defaults.stats.topoFDR  = 0;".
%                                                     NB: If it is left empty (i.e. [] or missing), then you will be asked to select it manually
%                               .p            <-- The threshold p-value. Should be a number, e.g. 0.05 or 0.001 or such.';
%                                                     NB:
%                                                        If it is left empty (i.e. [] or missing), then SPMthresStruct.MultTest_Cor will be checked
%                                                        if 'FWE'     then SPMthresStruct.p= 0.05;
%                                                        if 'none'    then SPMthresStruct.p= 0.001;
%                                                        if [](empty) then SPMthresStruct.p= []; (empty) and you will be asked to select it manually
%                               .k            <-- The minimum cluster size in voxels. Should be a number (integer>=0), e.g. 0 or 6 or 27 or 42 or 81 or such.';                 
%                                                     NB: If it is left empty (i.e. [] or missing), then you will be asked to select it manually
%  EffectiveDist_mm   <--   This indicates the MINIMUM distance in mm between local maxima. I.e. if two maxima are less than this distance apart,
%                           they will become one cluster with the higher one being the local maxima of the cluster.
%                           E.g. EffectiveDist_mm = 8;%mm OR EffectiveDist_mm = 16;%mm.
%                             NB: This does not mean that each cluster is maximally this size, clusters can become much larger
%                                 or even smaller based on the structure of the statistics map. 
%       
%  UseConnectMat      <--   This will activate the additional use of the connectedness matrix when creating clusters.
%                           This is a switch (boolean) variable, either 1 OR 0.
%                           The connectedness matrix indicates which voxels are connected via their neighbors,
%                           i.e. when "walking" from each voxel, which voxels can be reached when only using neighboring voxels,
%                           that are also significant and continuing from these to their neighbors and so on.
%                             NB1: If you want to keep SEPARATE "blobs" that are close to each other separated,
%                                  (i.e. "close" means closer than the "EffectiveDist_mm"),
%                                  then set UseConnectMat= 1; such that only connected regions can become singular cluster.
%                             NB2: This does not mean that connected areas will always be one singular cluster!!!
%                                  That still depends on the "EffectiveDist_mm".
%                             NB3: If you, to the contrary of NB1, want to combine separate "blobs" into one cluster,
%                                  i.e. if they are sufficiently close to each other, closer than the "EffectiveDist_mm",
%                                  then you should set UseConnectMat= 0; such that the clustering algorithm can combine them,
%                                  given the "EffectiveDist_mm".
%  ClusterNIIOutputDir <--  This (string) specifies the directory in which to save the clustering NIFTI-file.
%                           Clustering NIFTI-file will contain the the clusters indicated via integers from the most significant to the least significant.
%                             NB1: If this is left empty, then no file will be saved.
%                             NB2: If you want to create this file later, then just use the function LocMaxClusters2NIFTI.m,
%                                  see the help there on how to use it, using LocMaxStruct & MapExtractStruct.
%                             NB3: To avoid overwriting previous outputs, each file will have a date & time stamp in the filename.
%
%OUTPUTS:
%  LocMaxStruct         <--  The structure containing all information about the clustering. 
%                              NB: this will be cell with one entry per statistics values of negative and positive sign.
%  MapExtractStruct     <--  The structure containing all information about the extraction of significant statistics values
%                            from the input statistics image/contrast.
%  ClusterNIIOutputPath <-- The path to the clustering NIFTI-file that has been output. (if at all)
%
%
%USAGE:
%      [LocMaxStruct,MapExtractStruct,ClusterNIIOutputPath] = ClusterStatsMap(StatsMapPath,ThresholdSettings,EffectiveDist_mm,UseConnectMat,ClusterNIIOutputDir); 
%      [LocMaxStruct,MapExtractStruct] = ClusterStatsMap(StatsMapPath,ThresholdSettings,EffectiveDist_mm,UseConnectMat); %don't output clustering NIFTI
%      [LocMaxStruct,MapExtractStruct] = ClusterStatsMap(StatsMapPath,ThresholdSettings,EffectiveDist_mm); %don't use connections matrix %don't output clustering NIFTI
%
%
%
%V2.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment: V2.0(24.04.2015): Changed version of V1.5/6 ClusterStatsMap2Atlas.m; (02.February.2015): initial implementation based on test script.

%% %%%%%%%%%%%%%% setting for debug %%%%%%%%%%%%%%
%% Verbose?
% if(strcmp('Yes',questdlg('Output diagnostic messages?','Verbose?','Yes','No','Yes')))
%     Verbose = 1;
% else
%     Verbose = 0;
% end
Verbose = 0; %no diagnostics output

%% descriptive plots?
% if(strcmp('Yes',questdlg('Output descriptive/diagnostic plots of clustering iterations?','Descriptive plots?','Yes','No','Yes')))
%     DescripPlots = 1;
% else
%     DescripPlots = 0;
% end
DescripPlots = 0; %no descriptive plots for looking at iterations


%% check inputs
%% check inputs
%% check inputs

%% StatsMapPath: sensible input?
try
    if(isempty(StatsMapPath))
        StatsMapPath = spm_select(1,'any','Select statistics NIFTI OR SPM.mat for extracting statistics values for clustering...');
    end
    if(~ischar(StatsMapPath))
        if(iscellstr(StatsMapPath))
            StatsMapPath_tmp = StatsMapPath; clear StatsMapPath
            StatsMapPath = StatsMapPath_tmp{1};
        else
            error('"StatsMapPath" should be a string!');
        end
    end
catch CATCH_StatsMapPath
    disp_catch(CATCH_StatsMapPath,mfilename,'CATCH_StatsMapPath');
    StatsMapPath = spm_select(1,'any','Select statistics NIFTI OR SPM.mat for extracting statistics values for clustering...');
end
%% StatsMapPath: check the kind of input we have NIFTI or SPM.mat???
[StatsMapDir,StatsMapFName,StatsMapExt] = fileparts(StatsMapPath);
if(exist([StatsMapDir,filesep,StatsMapFName,StatsMapExt(1:4)],'file'))
    if(strcmpi(StatsMapFName,'SPM')&&strcmpi(StatsMapExt(1:4),'.mat'))
        StatsMap_isSPMmat = 1;
    else
        if(strcmpi(StatsMapExt(1:4),'.nii')||strcmpi(StatsMapExt(1:4),'.img')||strcmpi(StatsMapExt(1:4),'.hdr'))
            StatsMap_isSPMmat = 0;
        end
    end
else
    error(['"StatsMapPath"=="',StatsMapPath,'" is not a valid file!']);
end

%% ThresholdSettings: sensible input? & sensible in relation to "StatsMap_isSPMmat"?
try
    if(isempty(ThresholdSettings))
        if(StatsMap_isSPMmat)
            ThresholdSettings = struct('Contrast',[],'MultTestCor',[],'p',[],'k',[]);
            MapExtractStruct  = GetParamsFromSPMmat(StatsMapDir,ThresholdSettings); %contrast and threshold will be selected by user.
        else
            ThresholdSettings = nan(2,1);
            MapExtractStruct = GetParamsFromMap(StatsMapPath,ThresholdSettings); %in this case user will have to select thresholds
        end
    else
        if(StatsMap_isSPMmat)
            if(isstruct(ThresholdSettings))
                MapExtractStruct  = GetParamsFromSPMmat(StatsMapDir,ThresholdSettings); %contrast and threshold will either be selected automatically or if something is missing it will be selected by user.
            else
                disp('"ThresholdSettings" must be a structure if "StatsMapPath" points to a SPM.mat-file!');
                MapExtractStruct  = GetParamsFromSPMmat(StatsMapDir,[]); %contrast and threshold will be selected by user.
            end
        else
            MapExtractStruct = GetParamsFromMap(StatsMapPath,ThresholdSettings); %in case the thresholds are not defined user has to select them.
        end    
    end
catch CATCH_ThresholdSettings
    disp_catch(CATCH_ThresholdSettings,mfilename,'CATCH_ThresholdSettings');
    if(StatsMap_isSPMmat)
        ThresholdSettings = struct('Contrast',[],'MultTestCor',[],'p',[],'k',[]);
        MapExtractStruct  = GetParamsFromSPMmat(StatsMapDir,ThresholdSettings); %contrast and threshold will be selected by user.
    else
        ThresholdSettings = nan(2,1);
        MapExtractStruct = GetParamsFromMap(StatsMapPath,ThresholdSettings); %in this case user will have to select thresholds
    end
end

%% EffectiveDist_mm: sensible input?
try
    if(isempty(EffectiveDist_mm))
        %% minimum distance ???
        answer_SearchDist = inputdlg({'Effective/Search distance[mm]: '},'EffectiveDist_mm?',1,{'8'}); %get search distance from user
        EffectiveDist_mm  = eval(answer_SearchDist{1}); %search distance for local maxima
        clear answer_SearchDist %cleanup
    else
        if(~isnumeric(EffectiveDist_mm))
            answer_SearchDist = inputdlg({'Effective/Search distance[mm]: '},'EffectiveDist_mm?',1,{'8'}); %get search distance from user
            EffectiveDist_mm  = eval(answer_SearchDist{1}); %search distance for local maxima
            clear answer_SearchDist %cleanup
        else
            if(length(EffectiveDist_mm)>1)
                answer_SearchDist = inputdlg({'Effective/Search distance[mm]: '},'EffectiveDist_mm?',1,{'8'}); %get search distance from user
                EffectiveDist_mm  = eval(answer_SearchDist{1}); %search distance for local maxima
                clear answer_SearchDist %cleanup
            end
        end
    end
catch CATCH_EffectiveDist_mm
    disp_catch(CATCH_EffectiveDist_mm,mfilename,'CATCH_EffectiveDist_mm');
    answer_SearchDist = inputdlg({'Effective/Search distance[mm]: '},'EffectiveDist_mm?',1,{'8'}); %get search distance from user
    EffectiveDist_mm  = eval(answer_SearchDist{1}); %search distance for local maxima
    clear answer_SearchDist %cleanup
end

%% UseConnectMat: sensible input?
try
    if(isempty(UseConnectMat))
        if(strcmp('Yes',questdlg({'Use Connections matrix? (DEFAULT: "No")'; ' '; 'I.e. a matrix indicating for each voxel which other voxels (above threshold) are connected to it.'},'ConnectMat?','Yes','Def: No','Def: No')))
            UseConnectMat = 1;
        else
            UseConnectMat = 0;
        end
    else
        if(length(UseConnectMat)>1)
            if(strcmp('Yes',questdlg({'Use Connections matrix? (DEFAULT: "No")'; ' '; 'I.e. a matrix indicating for each voxel which other voxels (above threshold) are connected to it.'},'ConnectMat?','Yes','Def: No','Def: No')))
                UseConnectMat = 1;
            else
                UseConnectMat = 0;
            end
        end
    end
catch CATCH_UseConnectMat
    disp_catch(CATCH_UseConnectMat,mfilename,'CATCH_UseConnectMat');
    if(strcmp('Yes',questdlg({'Use Connections matrix? (DEFAULT: "No")'; ' '; 'I.e. a matrix indicating for each voxel which other voxels (above threshold) are connected to it.'},'ConnectMat?','Yes','Def: No','Def: No')))
        UseConnectMat = 1;
    else
        UseConnectMat = 0;
    end
end
        
%% ClusterNIIOutputDir: sensible input? --> if empty or not input --> don't use
try
    if(~isempty(ClusterNIIOutputDir))
        if(~exist(ClusterNIIOutputDir,'dir'))
            disp(['WARNING: Directory "',ClusterNIIOutputDir,'" does not exist.']);
            disp('WARNING: Will NOT write out Cluster NIFTI-file.');
            disp('NB: You can do this later manually using LocMaxClusters2NIFTI.m & the outputs of this function.');
            ClusterNIIOutputDir = [];
        end
    end
catch CATCH_ClusterNIIOutputDir
    disp_catch(CATCH_ClusterNIIOutputDir,mfilename,'CATCH_ClusterNIIOutputDir');
    disp('WARNING: Will NOT write out Cluster NIFTI-file.');
    disp('NB: You can do this later manually using LocMaxClusters2NIFTI.m & the outputs of this function.');
    ClusterNIIOutputDir = [];
end
if(~isempty(ClusterNIIOutputDir))
    if(StatsMap_isSPMmat)
        [tmp,fn] = fileparts(StatsMapDir); clear tmp
        OutFName= [fn,'_SPM_Con',num2str(MapExtractStruct.Positive.ConNr)];
    end
else
    OutFName = [];
end
    

    
%% %%%%%%%%%%     main program      %%%%%%%%%%
%% %%%%%%%%%%     main program      %%%%%%%%%%
%% %%%%%%%%%%     main program      %%%%%%%%%%

%% check which directions of data are available
DataPartsToDo = {};
ThresToDo     = [];
if(isfield(MapExtractStruct,'Negative'))
    if(~isempty(MapExtractStruct.Negative))
        DataPartsToDo{end+1,1} = 'Negative';
        ThresToDo(end+1) = MapExtractStruct.Thresholds(1);
    end
end
if(isfield(MapExtractStruct,'Positive'))
    if(~isempty(MapExtractStruct.Positive))
        DataPartsToDo{end+1,1} = 'Positive';
        ThresToDo(end+1,1) = MapExtractStruct.Thresholds(2);
    end
end
if(isempty(DataPartsToDo))
    error('No significant voxels to cluster!');
end

%% assign data (can be two collections of voxels for positive(1) and negative(2) stats-vals; NB: for negative "Local Maxima" is actually "Local-Minima")
LocMaxStruct = {}; %init empty
ClusterNIIOutputPath = {}; %init empty
for IndData = 1:length(DataPartsToDo)
    disp(['Treating "',DataPartsToDo{IndData},'" part of the input Map. (Thres>',num2str(ThresToDo(IndData)),').']);
    Coords    = MapExtractStruct.(DataPartsToDo{IndData}).Coords_mm; %dimensions are columns and datapoints are rows
    StatsVals = MapExtractStruct.(DataPartsToDo{IndData}).StatsVals; %statistic values per voxel
    VoxCoords = MapExtractStruct.(DataPartsToDo{IndData}).Coords_vox;
    
    %% use connectedness as well?
    if(UseConnectMat)
        %% Have a look at connectedness over neighbors
        tic
        ConnectionsMat = FindConnectedVoxels(VoxCoords,Verbose);
        t_FindConnectedVoxels = toc;
        disp(['Time needed to run "FindConnectedVoxels" is ',num2str(t_FindConnectedVoxels),'s']);
        ConnectionsSum = sum(ConnectionsMat,2);
        
        %% display connectedness
        if(DescripPlots)
            NVoxExp = 10^floor(log10((prod(quantile(ConnectionsSum(:),[.25 .5 .75])))^(1/length(quantile(ConnectionsSum(:),[.25 .5 .75]))))); %old; floor((median(ConnectionsSum(:))+quantile(ConnectionsSum(:),.25))/2); %median plus 20%
            figure(); clf;
            subplot(1,3,1); imagesc(ConnectionsMat,[0 1]); colormap(gray); title('connectivity matrix.');
            subplot(1,3,2); imagesc(repmat(ConnectionsSum,1,1),[0,NVoxExp]); title('sum of connections per voxel.');
            subplot(1,3,3); boxplot(ConnectionsSum,'notch','on','labels','ConnectionSums'); title('boxplot sum of connections.');
        end
        
        %% Move Voxels to local maxima
        LocMaxStruct_tmp = FindAllLocMax(Coords,StatsVals,EffectiveDist_mm,{'VoxelConnections',ConnectionsMat},{'Verbose',Verbose});
    else
        LocMaxStruct_tmp = FindAllLocMax(Coords,StatsVals,EffectiveDist_mm,{'Verbose',Verbose});
    end
    
    %% plots regarding iterations
    if(DescripPlots)
        %% fig for convergence plot
        if(evalin('base','exist(''Conv_fig_h'',''var'')'))
            Conv_fig_h = evalin('base','Conv_fig_h');
        else
            Conv_fig_h = figure();
            assignin('base','Conv_fig_h',Conv_fig_h);
        end
        
        %% show convergence
        figure(Conv_fig_h); clf;
        [AX,H1,H2]=plotyy(1:length(LocMaxStruct_tmp.Iterations.NClusters),LocMaxStruct_tmp.Iterations.NClusters,1:length(LocMaxStruct_tmp.Iterations.DistanceToConverge),LocMaxStruct_tmp.Iterations.DistanceToConverge);
        set(H1,'LineStyle','-','Marker','x');
        set(H2,'LineStyle','-','Marker','o');
        xlabel('Iteration'); title(['Convergence & number of Clusters for LocMax search [Minimum separation ',num2str(EffectiveDist_mm),'mm]']);
        % legend('Number of Clusters','Distance to CONVERGENCE');
        set(get(AX(1),'Ylabel'),'String','Number of Clusters')
        set(get(AX(2),'Ylabel'),'String','Distance to CONVERGENCE')
        
        %% change in distance matrix?
        if(evalin('base','exist(''ConvDist_fig_h'',''var'')'))
            ConvDist_fig_h = evalin('base','ConvDist_fig_h');
        else
            ConvDist_fig_h = figure();
            assignin('base','Conv_fig_h',ConvDist_fig_h);
        end
        
        Dist_ORG = squareform(pdist(Coords,'euclidean'),'tomatrix');
        
        figure(ConvDist_fig_h); clf;
        subplot(2,3,1); imagesc(Dist_ORG); title('Original distances'); %axis('square');
        for IndIter = 1:length(LocMaxStruct_tmp.Iterations.NClusters)
            Ind_i  = IndIter;
            Dist_i = LocMaxStruct_tmp.Iterations.DistMat_LocMax{IndIter};
            if((IndIter+1)<=length(LocMaxStruct_tmp.Iterations.NClusters))
                Ind_j = IndIter+1;
                Dist_j = LocMaxStruct_tmp.Iterations.DistMat_LocMax{IndIter+1};
            else
                Ind_j = IndIter;
                Dist_j = LocMaxStruct_tmp.Iterations.DistMat_LocMax{IndIter};
            end
            
            subplot(2,3,2); imagesc(Dist_i); title(['Distances  [    i= ',num2str(Ind_i),']']); %axis('square');
            subplot(2,3,3); imagesc(Dist_i-Dist_ORG,[-2.5*EffectiveDist_mm 2.5*EffectiveDist_mm]); title(['Difference [i-ORG=',num2str(Ind_j),'-ORG]']); %axis('square');
            
            subplot(2,3,4); imagesc(Dist_i); title(['Distances  [  i  =',num2str(Ind_i),']']); %axis('square');
            subplot(2,3,5); imagesc(Dist_j); title(['Distances  [j=i+1=',num2str(Ind_j),']']); %axis('square');
            subplot(2,3,6); imagesc(Dist_j-Dist_i,[-EffectiveDist_mm EffectiveDist_mm]); title(['Difference [ j-i =',num2str(Ind_j),'-',num2str(IndIter),']']); %axis('square');
            pause(1);
        end
    end
    
    %% make NIFTI!
    if(~isempty(ClusterNIIOutputDir))
        if(~StatsMap_isSPMmat)
            if(length(StatsMapExt)>4)
                OutFName = [StatsMapFName,'_Vol',StatsMapExt(6:end),'_',DataPartsToDo{IndData},'T',num2str(abs(ThresToDo(IndData)),3)];
            else
                OutFName = [StatsMapFName,'_',DataPartsToDo{IndData},'T',num2str(abs(ThresToDo(IndData)),3)];
            end
        end
        if(UseConnectMat)
            CMatStr= 'UseConnectMat';
        else
            CMatStr= '';
        end
        ClusterNIIOutputPath{end+1,1} = LocMaxClusters2NIFTI(LocMaxStruct_tmp.ClusterNo,MapExtractStruct.(DataPartsToDo{IndData}).Coords_vox,MapExtractStruct.V_map,[ClusterNIIOutputDir,filesep,'Cluster_',OutFName,'_',num2str(EffectiveDist_mm),'mmMinDist',CMatStr,'_',datestr(now,'yyyymmmdd_HHMMSS'),'.nii']);
        disp(' ');
    end
%     if(strcmp('Yes',questdlg('Create NIFTI from Clustering and then display Clustering as a check overview?','Display Clusters?','Yes','No','Yes')))
%         if(~exist('ConnectionsMat','var'))
%             [H,ColorsRGBperId,bg,OutputPath] = ShowOverlayClusters(LocMaxStruct_tmp.ClusterNo,VoxCoords,MapExtractStruct.V_map,[pwd,filesep,'Clusters_',num2str(EffectiveDist_mm),'mm_MinDist_',date]);
%         else
%             [H,ColorsRGBperId,bg,OutputPath] = ShowOverlayClusters(LocMaxStruct_tmp.ClusterNo,VoxCoords,MapExtractStruct.V_map,[pwd,filesep,'Clusters_',num2str(EffectiveDist_mm),'mm_MinDist_UsingConnectionsMatrix_',date]);
%         end
%     end
    
%% save LocMaxStruct?
%     if(strcmp('Yes',questdlg('Do you want to save the created LocalMaxima-Structure?','Save LocMaxStruct?','Yes','No','Yes')))
%         if(~SaveAll)
%             LocMaxStruct = rmfield(LocMaxStruct,'Iterations');
%         end
%         uisave({'LocMaxStruct'},['LocMaxStruct_FORstatsmapXYZ_',num2str(EffectiveDist_mm),'mm_MinDist.mat'])
%     end  
LocMaxStruct{end+1,1} = rmfield(LocMaxStruct_tmp,'Iterations'); %return current LocMaxStruct. BUT not each iteration, that would be too much.
    
end

%% Done
disp(' ');
disp('Done with ClusterStatsMap.');

end