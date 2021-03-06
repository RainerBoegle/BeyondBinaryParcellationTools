function [] = ComputeParcellationOfMaps(InputMapPaths,MapType,ThresholdType,Thresh,OutputDir)
% This function can be used to create parcellations of ICs/thresholded Maps/Normalized Cluster Maps.
% Just input the maps and indicate what type of input was used and specify a threshold, as well as an
% Output directory/Choice for generating an output directory in the input directory.
% 
% This function will create NIFTI-files containing information about the parcellation as well as
% MATLAB files containing the Parcellation Structure that contains other information and various
% matrices that can be used in further Graph Theoretical Analysis (GTA).
%
% List of NIFTI-files
% SumMap   --> All thresholded maps summed together, i.e. a value of "1" mean only one component at
%              this voxel and any other value means the corresponding number of overlaps.
% Base2Enc --> Base "2" Encoding of overlaps (NOT USED FOR DISPLAY, ONLY FOR LOOK UP)
%              Each component is associated with a power of 2 and summed together, -this results in
%              a unique assignment of overlaps.
% UniqueAssignment --> The Base2Encode map NIFTI is analysed and all unique overlaps are assign to a category. 
%                      This might result in a map that can be displayed with a colorbar for all clusters.
%                      (This is sometimes helpful)
% SingularICs   --> Only voxels that belong to one IC are assigned with the Number of the IC.
%                   (Useful for seeing which voxels primarily belong to one IC.)
% OverlapsOnly  --> The complement to the SingularICs NIFTI-file showing only those voxels where overlaps occur.
%                   (Sometimes helpful for seeing where particular overlaps occur.)
% WinnerTakeAll --> Show at each Voxel which IC has to highest value. 
%                   NB: THIS MIGHT ONLY YIELD REASONABLE RESULTS FOR "RAW" STATISTIC MAPS, NOT FOR THE "NORMALIZED CLUSTERS MAPS".
%
%
%
% ParcellationStruct
% .OverlapMat  -->   Matrix showing the overlap of ICs. The diagonal shows the number of voxels for
%                    each IC and all other entries show the number of overlaps between two ICs.
%                    (GTA of multiple entries can reveal multiple overlaps.)
% .AllLocMax   -->   Substructure containing all local maxima of all ICs and their co-occurrence
%  .ClusterNoAllICs (NVoxels,NICs)  --> Cluster Numbers 1 to NTotal many per IC (i.e. NTotal = max(ClusterNoAllICs(:)))
%  .ICnoPerCluster  (NTotal,1)      --> for each cluster list in which IC it is.
%  .NClusterPerIC   (NICs,1)        --> From IC 1 to IC N list the number of clusters. NTotal = sum(NClusterPerIC)
%  .LocMaxCoords_mm (NTotal,3)      --> Coordinate of the Local Maxima NTotal many
%  .OverlapMat      (NTotal,NTotal) --> Matrix showing the overlap of all clusters i.e. overlap of all possible clusters and therefore also the cluster size for each cluster.
%  .CoOccurrenceMat (NTotal,NTotal) --> Only those that are in the same IC are "co-occurring" here --> strucute of ICs will be visible
%  .CombCoOccurrenceMat (NTotal,NTotal) --> Combined Co-Occurrence Matrix i.e. if there is an overlap of the current cluster with another then mark it connected with the others that connect to the other cluster as well. 
%  NEED TO FIGURE OUT HOW TO REALLY DO THIS, i.e. to get all of this together and display it and so
%  on... (e.g. show the cluster matrix and then when interested in the clusters of a certain IC pick
%  all other clusters that show an overlap with these clusters or overlap through other ICs
%  clusters... damn difficult.) NEED TO LEARN GTA!!!

I need to change the way how clusters are treated, i.e. negative and positive numbers!!! Collected together using StatsSign

Also I need to make the input different for this function, i.e. I will need the original thresh_ maps and the results from clustering and normalization per cluster and so on...