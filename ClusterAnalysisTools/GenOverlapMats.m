function OverlapMats = GenOverlapMats(FilePaths1,FilePaths2,Thresh1,Thresh2)
% This function will generate overlap matrices for two sets of inputs (NIFTI-files)
% and thresholds for these inputs.
%
% The Overlap Matrices include
% 1.Dice Coefficient 
%               "DiceCoeff" = 2.*sum(MaskCl_1.*MaskCl_2)./(sum(MaskCl_1)+sum(MaskCl_2))
% 2.Minimum Cluster Overlap 
%               "MinClOverlap" = sum(MaskCl_1.*MaskCl_2)./min(sum(MaskCl_1),sum(MaskCl_2))
% 3.Overlap Size (diagonal is cluster size in voxels) 
%               "OverlapSize" = sum(MaskCl_1.*MaskCl_2)
% 4.Jaccard Coefficient (similar to dice coefficient but considers overlapping areas differently in numberator and denominator normalization.)
%               "JaccardCoeff" = sum(MaskCl_1.*MaskCl_2)./(sum(MaskCl_1)+sum(MaskCl_2)-sum(MaskCl_1.*MaskCl_2))
% 5.Sum Normalized Statistic Overlap in Cluster Mask, i.e., for a cluster all statistics values
%   are summed and each value divided by the sum, normalizing the sum to 1.
%   Therefore this overlap with another cluster, is the amount of statistics values that are
%   overlapping the other cluster while including the (normalized) significance amplitude! 
%   In contradistinction to the DiceCoeff which does not include significance amplitude!
%   NB: The other cluster is taken as a mask, i.e. its values are all set to 1!!! 
%       Not both taken as statistics values! 
%       --> Therefore this matrix is NOT SYMMETRIC!!!
%               "SumContribution" = sum(NormStatsValsCl_1.*MaskCl_2);
%   NB2: There are two versions of SumContribution12 indicating the normalized stats vals of FilePaths1 on the Masks created from FilePaths2  &  SumContribution21 indicating the normalized stats vals of FilePaths2 on the Masks created from FilePaths1.
% 6.Absolute value of the Correlation coefficient (Pearson correlation coefficient of statistics values!)
%               "AbsCorrCoeff" = abs(corr(StatsValsCl_1,StatsValsCl_2));
%
%
%Inputs:
%       FilePaths1,FilePaths2  (cellstring)  The paths to the statistics (or mask) files that should be used for OverlapMats creation.
%                                            If FilePaths1 contains N files and FilePaths2 contains M files then the Overlap Matrices will be N-x-M.
%       Thresh1,Thresh2        (Nfiles-x-2)  Negative and Positive thresholds for each input. If only one input is given then it will be used for  
%                                 Array      all files. If nothing is input or empty input is given, then the thresholds will be adjusted to take
%                                            all non-zero voxels.
%
%
%
%Usage:
%       OverlapMats = GenOverlapMats(FilePaths1,FilePaths2,Thresh1,Thresh2);
%       OverlapMats = GenOverlapMats(FilePaths1,FilePaths2); %take all nonzero values in the files indicated by Paths1 & Paths2
%       OverlapMats = GenOverlapMats(); %select files and take all nonzero values in the files
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (25.11.2016): initial implementation based on previous tests and GenerateOverlapMatrices.m

%% Check inputs
%FilePaths1
if(~exist('FilePaths1','var'))
    FilePaths1 = cellstr(spm_select(inf,'image','Select images for overlap matrix creation...'));
elseif(isempty(FilePaths1))
    FilePaths1 = cellstr(spm_select(inf,'image','Select images for overlap matrix creation...'));
else
    if(~iscellstr(FilePaths1))
        if(ischar(FilePaths1))
            FilePaths1 = cellstr(FilePaths1);
        else
            error('FilePaths1 must be a cellstring or convertible char-array!');
        end
    end
end
    
%FilePaths2
if(~exist('FilePaths2','var'))
    FilePaths2 = cellstr(spm_select(inf,'image','Select images for overlap matrix creation...'));
elseif(isempty(FilePaths2))
    FilePaths2 = cellstr(spm_select(inf,'image','Select images for overlap matrix creation...'));
else
    if(~iscellstr(FilePaths2))
        if(ischar(FilePaths2))
            FilePaths2 = cellstr(FilePaths2);
        else
            error('FilePaths2 must be a cellstring or convertible char-array!');
        end
    end
end
    
%Thresh1
if(~exist('Thresh1','var'))
    Thresh1 = repmat([-eps, eps],length(FilePaths1),1);
elseif(isempty(Thresh1))
    Thresh1 = repmat([-eps, eps],length(FilePaths1),1);
else
    if(isnumeric(Thresh1) && ((size(Thresh1,1)~=length(FilePaths1)) || (size(Thresh1,1)~=2)))
        if(size(Thresh1,1)~=length(FilePaths1))
            if((size(Thresh1,1)==1) && (size(Thresh1,2)==2))
                if((Thresh1(1)>0) && (Thresh1(2)<0))
                    Thresh1 = Thresh1([2, 1]);
                end
                Thresh1 = repmat(Thresh1,length(FilePaths1),1);
            elseif((size(Thresh1,1)==2) && (size(Thresh1,2)==1))
                Thresh1 = Thresh1';
                if((Thresh1(1)>0)&&(Thresh1(2)<0))
                    Thresh1 = Thresh1([2, 1]);
                end
                Thresh1 = repmat(Thresh1,length(FilePaths1),1);
            else
                error(['Thresh1 must be of size Nfiles(FilePaths1)==',num2str(length(FilePaths1)),'-x-2 OR 1-x-2(will be replicated)!']);
            end
        else
            if(size(Thresh1,1)~=2)
                error(['Thresh1 must be of size Nfiles(FilePaths1)==',num2str(length(FilePaths1)),'-x-2 OR 1-x-2(will be replicated)!']);
            end
        end
    else
        error(['Thresh1 must be NUMERIC array of size Nfiles(FilePaths1)==',num2str(length(FilePaths1)),'-x-2 OR 1-x-2(will be replicated)!']);
    end
end

%Thresh2
if(~exist('Thresh2','var'))
    Thresh2 = repmat([-eps, eps],length(FilePaths2),1);
elseif(isempty(Thresh1))
    Thresh2 = repmat([-eps, eps],length(FilePaths2),1);
else
    if(isnumeric(Thresh2) && ((size(Thresh2,1)~=length(FilePaths2)) || (size(Thresh2,1)~=2)))
        if(size(Thresh2,1)~=length(FilePaths2))
            if((size(Thresh2,1)==1) && (size(Thresh2,2)==2))
                if((Thresh2(1)>0) && (Thresh2(2)<0))
                    Thresh2 = Thresh2([2, 1]);
                end
                Thresh2 = repmat(Thresh2,length(FilePaths2),1);
            elseif((size(Thresh2,1)==2) && (size(Thresh2,2)==1))
                Thresh2 = Thresh2';
                if((Thresh2(1)>0)&&(Thresh2(2)<0))
                    Thresh2 = Thresh2([2, 1]);
                end
                Thresh2 = repmat(Thresh2,length(FilePaths2),1);
            else
                error(['Thresh2 must be of size Nfiles(FilePaths2)==',num2str(length(FilePaths2)),'-x-2 OR 1-x-2(will be replicated)!']);
            end
        else
            if(size(Thresh2,1)~=2)
                error(['Thresh2 must be of size Nfiles(FilePaths2)==',num2str(length(FilePaths2)),'-x-2 OR 1-x-2(will be replicated)!']);
            end
        end
    else
        error(['Thresh2 must be NUMERIC array of size Nfiles(FilePaths2)==',num2str(length(FilePaths2)),'-x-2 OR 1-x-2(will be replicated)!']);
    end
end

%% get the data as 4D Arrays X-x-Y-x-Z-x-NVols & X-x-Y-x-Z-x-MVols
disp('Loading data to create OverlapMats...');
[Data4D_1,Vols_1] = GetData4D(FilePaths1,Thresh1);
[Data4D_2,Vols_2] = GetData4D(FilePaths2,Thresh2);
N = length(Vols_1);
M = length(Vols_2);

%% create mask from both
disp('Creating mask from both inputs...');
Mask3D = double((sum(Data4D_1~=0,4)+sum(Data4D_2~=0,4))~=0);
NVoxels= length(find(Mask3D(:)~=0));

%% transform to 2D NVox(Mask)-x-NVols & NVox(Mask)-x-MVols
disp('Extracting data using Mask...');
Data2D_1 = zeros(NVoxels,N);
for Ind = 1:N
    CurrData = Data4D_1(:,:,:,Ind);
    Data2D_1(:,Ind) = CurrData(Mask3D(:)~=0);
end
Masks2D_1 = double(Data2D_1~=0);

Data2D_2 = zeros(NVoxels,M);
for Ind = 1:M
    CurrData = Data4D_2(:,:,:,Ind);
    Data2D_2(:,Ind) = CurrData(Mask3D(:)~=0);
end
Masks2D_2 = double(Data2D_2~=0);

%% Normalize StatsVals
disp('Normalize stats vals via total sum == 1...');
AbsData2D_1    = abs(Data2D_1);
SumAbsData2D_1 = sum(AbsData2D_1,1);
NormData2D_1   = AbsData2D_1./repmat(SumAbsData2D_1,NVoxels,1);

AbsData2D_2    = abs(Data2D_2);
SumAbsData2D_2 = sum(AbsData2D_2,1);
NormData2D_2   = AbsData2D_2./repmat(SumAbsData2D_2,NVoxels,1);

%% Generate Overlap Matrices via Matrix products
disp('Generating Overlap Matrices...');
ClusterSizes1   = repmat(sum(Masks2D_1,1)',1,size(Masks2D_2,2));
ClusterSizes2   = repmat(sum(Masks2D_2,1),size(Masks2D_1,2),1);

OverlapSize       =    Masks2D_1' * Masks2D_2;
DiceCoeff         = 2.*OverlapSize./(   ClusterSizes1+ClusterSizes2);
MinClOverlap      =    OverlapSize./min(ClusterSizes1,ClusterSizes2);
JaccardCoeff      =    OverlapSize./(   ClusterSizes1+ClusterSizes2-OverlapSize);
SumContribution12 =    NormData2D_1'* Masks2D_2;
SumContribution21 =    NormData2D_2'* Masks2D_1;
AbsCorrCoeff      = abs(corr(Data2D_1,Data2D_2));

%% assign struct
Inputs      = struct('FilePaths1',{FilePaths1},'FilePaths2',{FilePaths2},'Thresh1',Thresh1,'Thresh2',Thresh2,'N',N,'M',M,'NVoxels',NVoxels,'Mask3D',Mask3D,'Data2D_1',Data2D_1,'Data2D_2',Data2D_2);
OverlapMats = struct('Inputs',Inputs,...
                     'OverlapSize',OverlapSize,...
                     'DiceCoeff',DiceCoeff,...
                     'MinClOverlap',MinClOverlap,...
                     'JaccardCoeff',JaccardCoeff,...
                     'SumContribution12',SumContribution12,...
                     'SumContribution21',SumContribution21,...
                     'AbsCorrCoeff',AbsCorrCoeff);

%% Done.
disp(' ');
disp('Done.');
disp(' ');

end

%% subfunctions
function [Data4D,Vols] = GetData4D(FilePaths,Thresh)
% extract data

Vols = spm_vol(FilePaths);
Data4D = zeros(Vols{1}.dim(1),Vols{1}.dim(2),Vols{1}.dim(3),length(Vols));
for IndVol = 1:length(Vols)
    Data4D(:,:,:,IndVol) = (Vols{IndVol}.private.dat(:,:,:,Vols{IndVol}.n(1)).*double(Vols{IndVol}.private.dat(:,:,:,Vols{IndVol}.n(1))<=Thresh(IndVol,1))) + ...
                           (Vols{IndVol}.private.dat(:,:,:,Vols{IndVol}.n(1)).*double(Vols{IndVol}.private.dat(:,:,:,Vols{IndVol}.n(1))>=Thresh(IndVol,2)));
end

end