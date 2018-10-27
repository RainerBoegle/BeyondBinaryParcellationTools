function OverlapMatStruct = CreateOverlapMat(ThreshMapsList,ThreshType,Thresh) %OLD: [OverlapMat,ICnums,ThreshMapsList,NVoxelsTotal,Mask,Data2D,DataThresh2D,Thresh,DataDim] = CreateOverlapMat(ThreshMapsList,Thresh)
% This function will create the overlap matrix for all components given a list of NIFTI-files
% containing the path to the maps as a cell array (or string/char array).
%
% It will be assumed that the maps end in the IC number and this will be used to extract the
% numbers, if this fails then a warning will be given and you better know what you are doing and
% what the NIFTI-files are that you are inputing to get this matrix.
% I.e. if a numbering is found as expected then the inputs will be ordered according to it,
% if not then the inputs from the first to last will be used and numbering describes the input.
% (This is just a way of making data selection easier when dealing with the output of MELODIC stats.)
%
% NB: it is assumed that the maps are thresholded and ANY VOXEL NOT EQUAL TO ZERO WILL BE INCLUDED.
%
% The overlap matrix contains the number of shared voxels and can be displayed using the function
% "DispOverlapMatrix.m". 
% (NB2: the diagonal indicates the number of voxels for each respective input/IC.)
%
%Usage:
%      OverlapMatStruct = CreateOverlapMat(ThreshMapsList,ThreshType,Thresh);
%      OverlapMatStruct = CreateOverlapMat(ThreshMapsList); %ie take all.
%      OverlapMatStruct = CreateOverlapMat(ThreshMapsList,'P',0.2); %take only top 20%
%      OverlapMatStruct = CreateOverlapMat(ThreshMapsList,'Stats',[2.8 -Inf]); %take all above 2.8(a.u.) stats value and no negative ones (because -Inf)
%      OverlapMatStruct = CreateOverlapMat(ThreshMapsList,'Stats',[2.8   -3]); %take all above 2.8(a.u.) stats value and all below -3(a.u.) stats value.
%
%
%V1.5
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.5: (22.03.2016): allow different types of thresholding. V1.1: (06.03.2016): extension to collect data in structure for easy plotting V1.0: (20.02.2016): initial implementation

%% check inputs
%ThreshMapsList?
try
    if(~iscellstr(ThreshMapsList))
        if(ischar(ThreshMapsList))
            ThreshMapsList = cellstr(ThreshMapsList);
        end
    end
catch CATCH_ThreshMapsList
    disp_catch(CATCH_ThreshMapsList,'CreateOverlapMat>Input:ThreshMapsList','CATCH_ThreshMapsList');
    ThreshMapsList = cellstr(spm_select(Inf,'image','Select thresholded maps for creation of OverlapMat...'));
end

%ICnums?
[ok,ICnums] = findICnums(ThreshMapsList); %get IC nums assuming they are listed at the end of the file... (subfunction below)
if(ok) %could extract ICnums
    if(~issorted(ICnums)) %are they sorted?
        disp('IC# were found, BUT NOT SORTED! Will sort everything in ascending order...');
        [ICnums,sortInds] = sort(ICnums);
        ThreshMapsList = ThreshMapsList(sortInds); %ThreshMapsList must be a cellstr otherwise this doesn't work therefore the hassle above when checking the inputs.
    end
else
    disp('WARNING: could not determine ICnums, you better know what you are doing! ICnums will be left empty to reflect this fact.');
    ICnums = [];
end

%ThreshType & Threshold
if(~exist('ThreshType','var'))
    %Take all non-zero values! Use 'P' threshold because it is easier
    ThreshType = 'P';
    ThreshP    =  1; %use everything that is nonzero
    Thresh     =  1; %fix problems later
elseif(isempty(ThreshType))
    %Take all non-zero values! Use 'P' threshold because it is easier
    ThreshType = 'P';
    ThreshP    =  1; %use everything that is nonzero
    Thresh     =  1; %fix problems later
else
    switch(ThreshType)
        case 'P'
            ThreshP = Thresh;
        otherwise
            ThreshType  = 'Stats'; %just to be save.
            ThreshStats = Thresh;
    end
end    
%check & inform user
switch(ThreshType)
    case 'P'
        disp('Using threshold type "top X%" of voxels.');
        %check & inform user
        if(length(ThreshP)~=1)
            error('ThreshP must be a scalar!');
        elseif(ThreshP<=0)
            error('ThreshP must be >0 and <=1, i.e. not "top 0%" but between 100% & top X% (X>0)!');
        elseif(ThreshP>1)
            disp('ThreshP must be <=1 (& >0), will set it to 1.');
            ThreshP = 1;
        end
        if(ThreshP==1)
            disp('Using all significant voxels...');
        else
            disp(['Using the top ',num2str(ThreshP*100),'% of significant voxels for absolute maximum...']);
        end
    otherwise
        disp('Using threshold type statistics values (a.u.) above and below threshold.');
        ThreshType  = 'Stats'; %just to be save.
        if(isempty(Thresh))
            ThreshStats = [eps,-eps];
        else
            ThreshStats = Thresh;
        end
        %check & inform user
        if(length(ThreshStats)~=2)
            error('ThreshStats must two values!');
        else
            if(ThreshStats(1)<=0)
                error('ThreshStats(1) must be >0, i.e. positive threshold!');
            end
            if(ThreshStats(2)>=0)
                error('ThreshStats(2) must be <0, i.e. negative threshold!');
            end
        end
        if(all(isinf(ThreshStats))==1)
            disp('Using all significant voxels...');
        else
            if(ThreshStats(1))
                disp(['Using POSITIVE Voxels >',num2str(ThreshStats(1)),'(a.u.)...']);
            end
            if(ThreshStats(2))
                disp(['Using NEGATIVE Voxels <',num2str(ThreshStats(2)),'(a.u.)...']);
            end
        end
end

%% load data from NIFTI-files & Threshold
disp('Creating overlap matrix...');
Vtmp   = spm_vol(ThreshMapsList{1}); %assume equal size here
DataDim= Vtmp.dim;
Data2D = zeros(prod(Vtmp.dim),length(ThreshMapsList)); %Voxels-x-Inputs
Vols   = cell(length(ThreshMapsList),1); %all volume structures
%go over inputs and extract data
for Ind = 1:length(ThreshMapsList)
    Vtmp  = spm_vol(ThreshMapsList{Ind}); %assume equal size here
    if(DataDim(1)~=Vtmp.dim(1)||DataDim(2)~=Vtmp.dim(2)||DataDim(3)~=Vtmp.dim(3))
        error(['Data dimensions are not the same! (Check Input ',num2str(Ind),')']);
    else
        Vols{Ind} = Vtmp;
    end    
    DataTmp = squeeze(Vtmp.private.dat(:,:,:,Vtmp.n(1)));  
    switch(ThreshType)
        case 'P'
            Data2D(:,Ind)= DataTmp(:).*double((double(DataTmp(:)>(max(DataTmp(:))*(1-ThreshP)))+double(DataTmp(:)<(min(DataTmp(:))*(1-ThreshP))))>0); %threshold if needed.
        otherwise
            Data2D(:,Ind)= DataTmp(:).*double((double(DataTmp(:)>ThreshStats(1))+double(DataTmp(:)<ThreshStats(2)))>0); %threshold if needed.
    end
end
DataThresh2D = double(Data2D~=0);
Mask         = sum(Data2D~=0,2);
NVoxelsTotal = length(find(Mask~=0));

%% create matrix
OverlapMat = DataThresh2D'*DataThresh2D; %square matrix of #-Inputs x #-Inputs NB: the power of MATLAB! ;)

%% save in OverlapMatStruct
OverlapMatStruct.OverlapMat     = OverlapMat;
OverlapMatStruct.ThreshMapsList = ThreshMapsList;
OverlapMatStruct.ThreshType     = ThreshType;
OverlapMatStruct.Thresh         = Thresh;
OverlapMatStruct.ICnums         = ICnums;
OverlapMatStruct.NVoxelsTotal   = NVoxelsTotal;
%Data
OverlapMatStruct.Data.Vols        = Vols;
OverlapMatStruct.Data.DataDim     = DataDim;
OverlapMatStruct.Data.Mask        = Mask;
OverlapMatStruct.Data.DataThresh2D= DataThresh2D;
OverlapMatStruct.Data.Data2D      = Data2D;


%% Done.
disp(['Done with creation of OverlapMat(',num2str(size(OverlapMat,1)),',',num2str(size(OverlapMat,1)),')...']);

end