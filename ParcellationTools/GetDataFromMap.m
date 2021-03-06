function [XYZmm,StatsVals,XYZvox,Vol,VolInd] = GetDataFromMap(ThreshMapPath,Thresh)
% This function loads the data from the NIFTI file, 
% positive voxels are in the first cell and
% negative in the second entry of the cell.
%
%Usage:
%      [XYZmm,StatsVals,XYZvox,Vol] = GetDataFromMap(ThreshMapPath,Thresh);
%      [XYZmm,StatsVals,XYZvox,Vol] = GetDataFromMap(ThreshMapPath); %assume Thresh = [0 0];
%      [XYZmm,StatsVals,XYZvox,Vol] = GetDataFromMap(); %Select NIFTI and assume Thresh = [0 0];
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (16.03.2016): initial implementation based on previous tests

%% check inputs
%ThreshMapPath
if(~exist('ThreshMapPath','var'))
    ThreshMapPath = spm_select(1,'image','Select NIFTI file for extracting thresholded data...');
    disp(['Loading selected NIFTI-file "',ThreshMapPath,'"...']);
elseif(isempty(ThreshMapPath))
    ThreshMapPath = spm_select(1,'image','Select NIFTI file for extracting thresholded data...');
    disp(['Loading selected NIFTI-file "',ThreshMapPath,'"...']);
elseif(iscellstr(ThreshMapPath))
    ThreshMapPathTmp = ThreshMapPath{1}; clear ThreshMapPath
    ThreshMapPath = ThreshMapPathTmp; clear ThreshMapPathTmp
    disp(['Loading NIFTI-file: "',ThreshMapPath,'"...']);
elseif(ischar(ThreshMapPath))
    disp(['Loading "',ThreshMapPath,'"...']);
else
    error('"ThreshMapPath" must be a char/string or cellstring!');
end

%Thresh
if(~exist('Thresh','var'))
    disp('Assuming input map is thresholded, i.e. taking all voxels that are non-zero.');
    Thresh        = [0 0];
elseif(~isempty(Thresh))
    disp(['Will threshold map positively with Stats>',num2str(Thresh(1)),' and negatively with Stats<',num2str(Thresh(2)),'.']);
else
    disp('Assuming input map is thresholded, i.e. taking all voxels that are non-zero.');
    Thresh        = [0 0];
end


%% get volume struct
Vol = spm_vol(ThreshMapPath);

%% extract data 
Data3D = Vol.private.dat(:,:,:,Vol.n(1));
if(Vol.n(1)~=1)
    VolInd = Vol.n(1);
else
    VolInd = [];
end

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