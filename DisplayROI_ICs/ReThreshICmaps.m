function [ReThreshICmapFilePath,Vo_StatsFinal,ReThreshClustersFilePath,Vo_ClFinal] = ReThreshICmaps(ThreshICmapFilePath,ThreshPercentile,ClusterSizeThres)
% This function will use the Clustering Tools (taken from AtlasROITools V2.0) to cluster each input
% thresholded IC map (presumably thresholded via mixture model approach) and cluster these maps
% using a simple connectedness matrix clustering.
% Then each cluster is assumed as one dominant peak (with maybe some more subpeaks) and the
% statistic value differences ("stats-"distances" to peak") and distances from the peak are used
% to recluster each cluster such that each cluster gets more compact.
%
% This is done FOR EACH CLUSTER found in the NIFTI-file "ThreshICmapFilePaths" (string) in three (or four) steps:
% 1.by weighting the differences in statistics values exponentially while the distance to
%   the peak is weighted linearly
% 2.both weightings are multiplied.
% 3.The resulting distance values (product of weighting of "stats-distances" and voxel-distances)
%   are then thresholded using the input Percentile, i.e. throwing out those voxels in the cluster
%   that are "further away" in terms of the final distance than the Percentile for thresholding.
%   E.g. if the Percentile for thresholding (Input "ThreshPercentile") is 0.5 (50%) then the
%   threshold is effectively the median (median is the 50%ile).
%   THE DEFAULT VALUE FOR "ThreshPercentile" IS 0.5, I.E. USE THE MEDIAN.
% 4.If the remaining cluster size is then below the cluster threshold (Input "ClusterSizeThres"), then it will be removed.
%   THE DEFAULT VALUE FOR "ClusterSizeThres" IS 13, I.E. HALF OF A FULL NEIGHBORHOOD OF 1. 
%   NB: A neighborhood of 1 means a voxel with all its nearest neighbors, i.e. 27. 
%       Therefore, half the full neighborhood is 13==floor(27/2).
%
%
%Usage:
%       [ReThreshICmapFilePath] = ReThreshICmaps(ThreshICmapFilePath,ThreshPercentile,ClusterSizeThres);       
%       [ReThreshICmapFilePath] = ReThreshICmaps(ThreshICmapFilePath,0.5,13); %Rethreshold the clusters found in "ThreshICmapFilePaths" (string), using median of distance values as threshold above which voxels are excluded and a final cluster size of 13 voxels, see above.
%       [ReThreshICmapFilePath] = ReThreshICmaps(ThreshICmapFilePath);        %Same as before, but "ThreshPercentile" is set to 0.5 and "ClusterSizeThres" is set to 13, BY DEFAULT.
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (08.11.2015): initial implementation (based on test script version "ATestInThresholdingAfterClustering.m")

%% Settings
PlotData = 0; %plot results & intermediate steps

%% Check inputs
%ThreshICmapFilePath
try
    if(isempty(ThreshICmapFilePath))
        ThreshICmapFilePath = spm_select(1,'image','Select Thresholded IC map that should be rethresholded...');
    else
        if(~ischar(ThreshICmapFilePath))
            error('"ThreshICmapFilePath" has to be a string!');
        end
    end
catch CATCH_ThreshICmapFilePath
    disp_catch(CATCH_ThreshICmapFilePath,[mfilename,'>ThreshICmapFilePath'],'CATCH_ThreshICmapFilePath');
    ThreshICmapFilePath = spm_select(1,'image','Select Thresholded IC map that should be rethresholded...');
end
            
%ThreshPercentile
try
    if(isempty(ThreshPercentile))
        disp('"ThreshPercentile" is not set! Will default to ThreshPercentile=0.5.');
        ThreshPercentile = 0.5;
    else
        if((ThreshPercentile<=0)||(ThreshPercentile>=1))
            disp(['"ThreshPercentile" must be between 0+ and 1-, i.e. ThreshPercentile must be from open interval (0,1)! [Input=',num2str(ThreshPercentile),']']);
        end
    end
catch CATCH_ThreshPercentile
    disp_catch(CATCH_ThreshPercentile,[mfilename,'>ThreshPercentile'],'CATCH_ThreshPercentile');
    disp('"ThreshPercentile" is not set! Will default to ThreshPercentile=0.5.');
    ThreshPercentile = 0.5;
end

%ClusterSizeThres
try
    if(isempty(ClusterSizeThres))
        disp('"ClusterSizeThres" is not set! Will default to ClusterSizeThres=13.');
        ClusterSizeThres = 13;
    else
        if(ClusterSizeThres==0);
            disp('"ClusterSizeThres" is set to zero. Will skip exclusion of clusters by their size. I hope you better know what are doing...');
        else
            if(ClusterSizeThres<0)
                error(['"ClusterSizeThres" must be greater or equal zero! [Input=',num2str(ClusterSizeThres),']']);
            end
        end
    end
catch CATCH_ClusterSizeThres
    disp_catch(CATCH_ClusterSizeThres,[mfilename,'>ClusterSizeThres'],'CATCH_ClusterSizeThres');
    disp('"ClusterSizeThres" is not set! Will default to ClusterSizeThres=13.');
    ClusterSizeThres = 13;
end

%% get ThresIC Volume info
disp(' ');
disp(['Loading "',ThreshICmapFilePath,'" for rethresholding clusters...']);

VthresIC = spm_vol(ThreshICmapFilePath);
StatsData= VthresIC.private.dat(:);
StatsData(StatsData==0) = []; %remove non-significant

Thresholds    = nan(2,1);
if(any(StatsData<0))
    Thresholds(1) = max(StatsData(StatsData<0));
end
if(any(StatsData>0))
    Thresholds(2) = min(StatsData(StatsData>0));
end
if(all(isnan(Thresholds)))
    error('No nonzero Voxels found!');
else
    disp(['Will use thresholds [',num2str(Thresholds'),'] for "',ThreshICmapFilePath,'".']);
end
 
%% load data from selected IC
[ThresIC_MapExtractStruct] = GetParamsFromMap(ThreshICmapFilePath,Thresholds);

%% find connect mat & clusters from stats vals
ThresDirs    = {'Positive';'Negative'};
ConnectMat   = cell(length(ThresDirs),1);
ClusterNr_vox= cell(length(ThresDirs),1);
for IndThresDirection = 1:length(ThresDirs)
    if(~isempty(ThresIC_MapExtractStruct.(ThresDirs{IndThresDirection})))
        disp(ThresDirs{IndThresDirection})
        %% make connectedness matrix
        ConnectMat{IndThresDirection} = FindConnectedVoxels(ThresIC_MapExtractStruct.(ThresDirs{IndThresDirection}).Coords_vox);
        
        %% find clusters within this ranked by significance of peak
        ClusterNr_vox{IndThresDirection} = ConnectMat2Clusters(ConnectMat{IndThresDirection},ThresIC_MapExtractStruct.(ThresDirs{IndThresDirection}).StatsVals);
    end
end

%% collect all data together First Positive then Negative.(if available)
ClusterNr_voxTotal = [];
StatsVals_Total    = [];
Coords_voxTotal    = [];
for IndThresDirection = 1:length(ThresDirs)
    if(~isempty(ThresIC_MapExtractStruct.(ThresDirs{IndThresDirection})))
        disp(ThresDirs{IndThresDirection})
        Shift = max(ClusterNr_voxTotal); 
        if(isempty(Shift)) 
            Shift = 0; 
        end
        ClusterNr_voxTotal = [ClusterNr_voxTotal; ClusterNr_vox{IndThresDirection}+Shift];
        StatsVals_Total    = [StatsVals_Total;    ThresIC_MapExtractStruct.(ThresDirs{IndThresDirection}).StatsVals];
        Coords_voxTotal    = [Coords_voxTotal;    ThresIC_MapExtractStruct.(ThresDirs{IndThresDirection}).Coords_vox];
    end
end

%% rethreshold clusters
StatsVals_Final   = StatsVals_Total;
UniqueClusterNums = unique(ClusterNr_voxTotal);
for IndCl = 1:length(UniqueClusterNums)
    disp(['Cluster ',num2str(UniqueClusterNums(IndCl))]);
    %prep indices and data for rethresholding
    IDXoI        = find(ClusterNr_voxTotal==UniqueClusterNums(IndCl));
    StatsVals_oI = StatsVals_Total(IDXoI);
    Coords_oI    = Coords_voxTotal(IDXoI,:);
    if(all(StatsVals_oI>=0))
            [PeakStats, IdxPeak] = max(StatsVals_oI(:));
        elseif(all(StatsVals_oI<=0))
            [PeakStats, IdxPeak] = min(StatsVals_oI(:));
        else
            error(['Mixed stats for Cluster ',num2str(UniqueClusterNums(IndCl))]);
    end
    PeakCoord = Coords_oI(IdxPeak,:);
        
    %calculate & reweight distance and stats-"distance"
    VoxDist   = sqrt(sum((Coords_oI-repmat(PeakCoord,size(Coords_oI,1),1)).^2,2)); VoxDist = 100.*VoxDist./max(VoxDist);
    StatsDist = abs(StatsVals_oI-PeakStats);                                     StatsDist = 100.*StatsDist./max(StatsDist);
    
    %dVal = (StatsDist.^3).*VoxDist; %quadratically grow StatsDistance to make it more important than the distance, -which enters linearly.
    dVal = (exp(StatsDist)-1).*VoxDist;  %exponentially grow StatsDistance to make it more important than the distance, -which enters linearly.
    %ThreshPercentile = 0.5; %1/exp(1); %median for 0.5 otherwise X==%percentile NB: 1/exp(1) is ~37%ile
    CurrThresh = quantile(dVal,ThreshPercentile); %current threshold (for this cluster) is the percentile defined above and we only take those below it.
    StatsVals_Final(IDXoI(dVal>CurrThresh)) = 0; %set those that are "too far away" to zero.
    
    if(PlotData)
        figure();
        subplot(2,2,[1 2]); plot(VoxDist( dVal<CurrThresh),StatsDist( dVal<CurrThresh),'kx'); title(['Cluster ',num2str(UniqueClusterNums(IndCl))]); xlabel('VoxDist'); ylabel('StatsDist'); hold on
        subplot(2,2,[1 2]); plot(VoxDist(dVal>=CurrThresh),StatsDist(dVal>=CurrThresh),'r*'); legend('retained','removed');
        subplot(2,2,   3);  plot(VoxDist,  dVal./max(dVal),'kx'); xlabel('VoxDist');   ylabel('dVal'); hold on;
        subplot(2,2,   3);  plot(  0:100,    ThreshPercentile,'b--');
        subplot(2,2,   4);  plot(StatsDist,dVal./max(dVal),'kx'); xlabel('StatsDist'); ylabel('dVal'); hold on
        subplot(2,2,   4);  plot(    0:100,  ThreshPercentile,'b--');
    end
end
ClusterNr_voxFinal = ClusterNr_voxTotal.*(StatsVals_Final~=0); %remove cluster.

%% evaluate cluster sizes
ClusterSizeThres = 13; %half of full neighborhood of 1, i.e. floor(27/2). %7; %simple neighborhood of 1. %27; %full neighborhood of 1. %0; %threshold for cluster size if 0 then no threshold is applied
if(ClusterSizeThres)
    UniqueClusterNumsFinal = unique(ClusterNr_voxFinal);
    SizeUniqueClusterNumsFinal = zeros(length(UniqueClusterNumsFinal),1);
    for Ind = 1:length(UniqueClusterNumsFinal)
        SizeUniqueClusterNumsFinal(Ind) = length(find(ClusterNr_voxFinal==UniqueClusterNumsFinal(Ind)));
        if(SizeUniqueClusterNumsFinal(Ind)<ClusterSizeThres)
            ClusterNr_voxFinal(ClusterNr_voxFinal==UniqueClusterNumsFinal(Ind)) = 0;
        end
    end
    StatsVals_Final(ClusterNr_voxFinal==0) = 0;
    ClusterSizeThresStr = ['CLsizeThres',num2str(ClusterSizeThres)];
    
    disp([num2str(max(UniqueClusterNumsFinal)),' Clusters have been assigned. (ClusterSizes are [',num2str(SizeUniqueClusterNumsFinal'),'] Voxels.)']);
    disp(['Eliminating those below ',num2str(ClusterSizeThres),' Voxels.']);
else
    disp([num2str(max(UniqueClusterNumsFinal)),' Clusters have been assigned. (ClusterSizes are [',num2str(SizeUniqueClusterNumsFinal'),'] Voxels.)']);
    ClusterSizeThresStr = '';
end

%% output & display(if wanted
[OutDir,Fname,ext] = fileparts(ThreshICmapFilePath);
if(isempty(OutDir))
    OutDir = pwd;
end

%rethresholded stats.
[ReThreshICmapFilePath,Vo_StatsFinal] = LocMaxClusters2NIFTI(StatsVals_Final,Coords_voxTotal,VthresIC,[OutDir,filesep,'ReThresh_',Fname,ext]);

%rethresholded clusters.
[ReThreshClustersFilePath,Vo_ClFinal] = LocMaxClusters2NIFTI(ClusterNr_voxFinal,Coords_voxTotal,VthresIC,[OutDir,filesep,'ReThreshClusters_',Fname,ext]);

if(PlotData)
    %original
    [H,SlObj,params,SliceIndices] = DisplayStats('Stats',ThreshICmapFilePath); clear H SlObj params
    H = helpdlg(['Displaying ORIGINAL stats-map "',ThreshICmapFilePath,'"'],['"',ThreshICmapFilePath,'"']);
    uiwait(H);
    
    %rethresholded stats.
    DisplayStats(SliceIndices,ReThreshICmapFilePath);
    H = helpdlg(['Displaying RETHRESHOLDED stats-map saved in file "ReThresh_',Fname,ext],['"ReThresh_',Fname,ext]);
    uiwait(H);
    
    %rethresholded clusters.
    DisplayClusters(SliceIndices,ReThreshClustersFilePath);
    H = helpdlg(['Displaying ReThreshClusters: "',ReThreshClustersFilePath,'"'],['"',ReThreshClustersFilePath,'"']);
    uiwait(H);
end

%% Done.
disp('...done');
disp(' ');

end
