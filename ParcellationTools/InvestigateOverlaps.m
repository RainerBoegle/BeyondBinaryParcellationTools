function [FilesOI,ICsOI,VoxOI] = InvestigateOverlaps(CompEncodeDir,CoordsOI,Nhood,LimitSumMap)
% This function can be used to extract the overlap information using the component encoding maps.
% Specifically this concerns the "SumMap.nii" & "Base2Encode.nii" Maps.
% 
% Given coordinates in "CoordsOI" in [mm] MNI-space and the size of the Neighborhood "Nhood" (number of
% neighboring voxels to include) and the Limit for SumMap (i.e. how many overlaps have to
% occur in a voxel such that it is included in the analysis), this function will read the overlaps
% from "SumMap.nii" (for the coordinate & Neighborhood) and selects those voxels that are above
% the Limit for the SumMap i.e. the number of overlaps.
% These voxels that were selected in this first step from "SumMap.nii" are then used for reading
% the values from the Map "Base2Encode.nii" that allows an unambiguous identification of the ICs
% that mix together in each of the voxels.
%
% All unique ICs for a certain Coordinate are collected and will be returned.
% With these it is possible to select the input data paths from "DataStruct"
% FilesOI = DataStruct.InputFiles(DataStruct.ICnums==ICsOI);
% These can then be plotted.
%
%FUTURE EXTENSION:
% 1. find all clusters (connected voxels) in "SumMap.nii" with limit>=LimitSumMap
% 2. for each cluster find all ICsOI
% 3. get all the clusters IN ALL the selected ICsOI (maybe name according to atlas?)
% 4. make a co-occurrence matrix over all clusters in all the ICsOI.
%    (of course these link over the overlap areas that we started with, but it is intersting to know
%     which ones these are and maybe we can automatically ascribe them.)
% 3b&4b.maybe try to see which clusters across the ICsOI do overlap such that they can be merged???
%       one measure would be overlap of voxels and if LocalMaxima of clusters are contained in the
%       extend of the other cluster, i.e. clusters do not need to overlap perfectly, but their local
%       maxima might be shared voxels of both???
%
%
%Inputs: 
%       CompEncodeDir (string)     Path to the directory that contains the "SumMap.nii", "Base2Encode.nii" & "DataStruct.mat" files
%       CoordsOI      (NVoxOI-x-3) The coordinates to be investigated.
%       Nhood         (1-x-1)      The size of the Neighborhood in Voxels that should be covered.
%       LimitSumMap   (1-x-1)      The limit i.e. lowest possible value in "SumMap.nii" == i.e. number of overlaps 
%                                  that a voxel should have to be included in the analysis of "Base2Encode.nii",
%                                  for identifying the ICs of interest (ICsOI).
%
%
%Usage:
%      [FilesOI,ICsOI,VoxOI] = InvestigateOverlaps(CompEncodeDir,CoordsOI,Nhood,LimitSumMap);
%      [FilesOI,ICsOI,VoxOI] = InvestigateOverlaps(CompEncodeDir,[0 -30 40; 0 -29 50;6 -30 50; -8 -30 49],2,3);
%
%
%
%V0.5
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V0.5: (22.03.2016): initial implementation WITHOUT SAFETY CHECKS FOR INPUTS

%% Check Inputs

%% Get data
V_SumMap   = spm_vol([CompEncodeDir,filesep,'SumMap.nii']);
V_Base2Enc = spm_vol([CompEncodeDir,filesep,'Base2Encode.nii']);
load([CompEncodeDir,filesep,'DataStruct.mat']);

%% do tranfo of MM-coords(MNI) to Voxel-Coords
mm2v = inv(V_SumMap.mat);
TestDiff = mm2v(1:3,1:3)-diag(diag(mm2v(1:3,1:3)));
if(all(TestDiff(:)==0)) %rotation&scaling submatrix is diagonal --> we can do it fast with bsxfun
    CoordsOI_Vox=round(bsxfun(@plus,bsxfun(@times,CoordsOI,diag(mm2v(1:3,1:3))'),mm2v(1:3,4)'));
else %not diagonal
    CoordsOI_Vox=zeros(size(CoordsOI,1),3);
    for k = 1:size(CoordsOI,1)
        CoordsOI_Vox(k,:) = CoordsOI(k,:)*mm2v(1:3,1:3)+mm2v(1:3,4)';
    end
    CoordsOI_Vox = round(CoordsOI_Vox);
end

%% go over coords and get information
FilesOI = cell(size(CoordsOI,1),1);
ICsOI   = cell(size(CoordsOI,1),1);
VoxOI   = cell(size(CoordsOI,1),1);

disp('Collecting all coords and data for finding ICsOI and saving the paths...');
for IndCoord = 1:length(FilesOI)
    %% create possible vox-coords
    CurrCoord = CoordsOI_Vox(IndCoord,:);
    AllPossibleCoords = GetNhood(CurrCoord,V_SumMap.dim,Nhood);
    
    %% get data (SumMap) from coords
    Data3D = V_SumMap.private.dat(:,:,:);
    LinInds_AllPossibleCoords = sub2ind(V_SumMap.dim,AllPossibleCoords(:,1),AllPossibleCoords(:,2),AllPossibleCoords(:,3));
    SelectedData = Data3D(LinInds_AllPossibleCoords);
    
    %% filter out those values below LimitSumMap
    AllPossibleCoords(SelectedData<LimitSumMap,:) = []; %remove
    LinInds_AllPossibleCoords = sub2ind(size(Data3D),AllPossibleCoords(:,1),AllPossibleCoords(:,2),AllPossibleCoords(:,3));
    VoxOI{IndCoord} = AllPossibleCoords; %save coords that are according to all requirements
    
    %% get data from Base2Encode & find ICsOI
    Data3D = V_Base2Enc.private.dat(:,:,:);
    SelectedData = Data3D(LinInds_AllPossibleCoords);
    
    ICcollection = [];
    for Ind = 1:length(SelectedData)
        tempBinary = arrayfun(@str2num,dec2bin(SelectedData(Ind))); %take current component mix number (base2encode) make binary and convert to a vector of 0&1
        tempBinary = tempBinary(end:-1:1); %reverse such that order is right
        CurrICsOI  = find(tempBinary~=0);
        ICcollection = [ICcollection; CurrICsOI(:)]; %find those positions with a "1" --> these are the IC numbers!
    end
    ICsOI{IndCoord} = unique(ICcollection); %save unique ICs in the collection
    
    %% get files from DataStruct
    AllFilesCurr = cell(length(ICsOI{IndCoord}),1);
    for Ind = 1:length(ICsOI{IndCoord})
        AllFilesCurr{Ind} = DataStruct.InputFiles{DataStruct.ICnums==ICsOI{IndCoord}(Ind)};
    end
    FilesOI{IndCoord} = AllFilesCurr;
end

%% Done.
disp('DONE.');
disp(' ');

end

%% subfunctions
%% GetNhood
function AllPossibleCoords = GetNhood(CurrCoord,Dim,Nhood)
%collect all possible voxel coords

AllPossibleCoords      = zeros((2*Nhood+1)^3,3);
AllPossibleCoords(1,:) = CurrCoord;
Steps = (-Nhood:Nhood);

Index = 0; %init
IndsRemove = []; %init empty
for IndZ = Steps+CurrCoord(3)
    for IndY = Steps+CurrCoord(2)
        for IndX = Steps+CurrCoord(1)
            Index = Index+1;
            if((IndX<=Dim(1)&&IndX>=1)&&(IndY<=Dim(2)&&IndY>=1)&&(IndZ<=Dim(3)&&IndZ>=1))
                if((IndX==CurrCoord(1))&&(IndY==CurrCoord(2))&&(IndZ==CurrCoord(3))) %remove this one because we moved it to the beginning.
                    IndsRemove = [IndsRemove; Index];
                else
                    AllPossibleCoords(Index,:) = [IndX IndY IndZ];
                end
            else
                IndsRemove = [IndsRemove; Index];
            end
        end
    end
end

%% remove those that were outside if necessary
if(~isempty(IndsRemove))
    AllPossibleCoords(IndsRemove,:) = []; %remove
    AllPossibleCoords = [CurrCoord; AllPossibleCoords];
end

end
