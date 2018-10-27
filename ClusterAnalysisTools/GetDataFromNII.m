function [StatsVals,XYZvox,XYZmm] = GetDataFromNII(Vol,Thresh)
% This function loads the data from the NIFTI file, 
% FIRST the positive voxels and THEN the negative voxels, that SURVIVE thresholding.
% 
% If a 4D-file is selected in its entirety then the data is returned as
% cell-vectors each corresponding to each of the volumes.
%
%Usage:
%       [StatsVals,XYZvox,XYZmm] = GetDataFromNII(Vol,Thresh);
%       [StatsVals,XYZvox,XYZmm] = GetDataFromNII(Vol); %get everything non-zero.
%       [StatsVals,XYZvox,XYZmm] = GetDataFromNII(spm_vol(spm_select(inf,'image','Select NIFTI-file...'))); %Select file using spm_select and get volume via spm_vol.
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (01.11.2016): initial implementation based on previous tests

%% check inputs
%Thresh
if(~exist('Thresh','var'))
    Thresh = [eps, -eps];
else
    if(isempty(Thresh))
        Thresh = [eps, -eps];
    else
        if((sign(Thresh(1))==-1) && ((sign(Thresh(2))==+1))) %need to flip it
            Thresh = Thresh(2:-1:1);
        else
            if(((sign(Thresh(1))==-1) && (sign(Thresh(2))==-1)) || ((sign(Thresh(1))==+1) && (sign(Thresh(2))==+1)))
                error('Thresh should have a positive and a negative value or empty if none used.');
            end                
        end
    end
end

%Vol 
if(~exist('Vol','var'))
    Vol = spm_vol(spm_select(inf,'image','Select NIFTI-file(s)...'));
end
if(length(Vol)>1) %USE RECURSIVE CALLS!
    StatsVals  = cell(length(Vol),1);
    XYZvox     = cell(length(Vol),1);
    XYZmm      = cell(length(Vol),1);
    reverseStr = ''; %init
    if(isstruct(Vol))
        for IndVol = 1:length(Vol)
            msg = sprintf(['Extracting Volume %0',num2str(ceil(log10(length(Vol)))),'d of %0',num2str(ceil(log10(length(Vol)))),'d...'], IndVol, length(Vol));
            fprintf([reverseStr, msg]); %delete the last message and print current message.
            reverseStr = repmat(sprintf('\b'), 1, length(msg)); %setup the next reverse string
            if(nargout>2) %Calculate only if mm-Coords are needed.
                [StatsVals{IndVol},XYZvox{IndVol},XYZmm{IndVol}] = GetDataFromNII(Vol(IndVol),Thresh);
            else
                [StatsVals{IndVol},XYZvox{IndVol}] = GetDataFromNII(Vol(IndVol),Thresh);
            end
        end
        msg = sprintf(['Extracting Volume %0',num2str(ceil(log10(length(Vol)))),'d of %0',num2str(ceil(log10(length(Vol)))),'d done.'], IndVol, length(Vol));
        fprintf([reverseStr, msg]); %delete the last message and print current message.
        disp(' ');
    else
        if(iscell(Vol))
            for IndVol = 1:length(Vol)
                msg = sprintf(['Extracting Volume %0',num2str(ceil(log10(length(Vol)))),'d of %0',num2str(ceil(log10(length(Vol)))),'d...'], IndVol, length(Vol));
                fprintf([reverseStr, msg]); %delete the last message and print current message.
                reverseStr = repmat(sprintf('\b'), 1, length(msg)); %setup the next reverse string
                if(nargout>2) %Calculate only if mm-Coords are needed.
                    [StatsVals{IndVol},XYZvox{IndVol},XYZmm{IndVol}] = GetDataFromNII(Vol{IndVol},Thresh);
                else
                    [StatsVals{IndVol},XYZvox{IndVol}] = GetDataFromNII(Vol{IndVol},Thresh);
                end
            end
            msg = sprintf(['Extracting Volume %0',num2str(ceil(log10(length(Vol)))),'d of %0',num2str(ceil(log10(length(Vol)))),'d done.'], IndVol, length(Vol));
            fprintf([reverseStr, msg]); %delete the last message and print current message.
            disp(' ');
        end
    end
    return;
else
    if(isempty(Vol))
        error('Empty volume input!');
    else
        if(iscell(Vol)) %careful this could be due to spm_vol's odd behavior when a path to a 4D-NIFTI is input as one cell --> check contents
            if(length(Vol{1})~=1||isstruct(Vol{1}))
                Vol = Vol{1};
                [StatsVals,XYZvox,XYZmm] = GetDataFromNII(Vol,Thresh); %recursive call!
                return;
            end
        end
    end
end

%% The actual processing...
%% extract data 
Data3D = Vol.private.dat(:,:,:,Vol.n(1));

%% assign data accordingly
XYZmmTmp     = cell(2,1);
StatsValsTmp = cell(2,1);
XYZvoxTmp    = cell(2,1);

%positive stats vals
if(any(Data3D(:)>Thresh(1)))
    StatsValsTmp{1} = Data3D(Data3D(:)>Thresh(1)); %statistics values POSITIVE
    IndsStatsVals   = find(Data3D(:)>Thresh(1));   %corresponding indices
    XYZvoxTmp{1}    = zeros(length(IndsStatsVals),3); %nVoxels-x-3 (spatial dimensions)
    [XYZvoxTmp{1}(:,1),XYZvoxTmp{1}(:,2),XYZvoxTmp{1}(:,3)] = ind2sub(Vol.dim,IndsStatsVals); %indices to 3D/subscript
    
    if(nargout>2) %Calculate only if mm-Coords are needed.
        %do trafo from voxels to world (mm-coords)
        TestDiff = Vol.mat(1:3,1:3)-diag(diag(Vol.mat(1:3,1:3))); %if this contains only zeros then we can do this very fast!
        if(all(TestDiff(:)==0)) %rotation&scaling submatrix is diagonal --> we can do it fast with bsxfun
            XYZmmTmp{1}=bsxfun(@plus,bsxfun(@times,XYZvoxTmp{1},diag(Vol.mat(1:3,1:3))'),Vol.mat(1:3,4)');
        else %not diagonal
            disp('WARNING(pos): VOXEL2WORLD 3D submatrix is NOT diagonal!!! (NB: usually any normalized (MNI) space should be diagonal in the 3D subspace.)'); %usually any normalized (MNI) space should be diagonal in the 3D subspace.
            XYZmmTmp{1}=zeros(length(IndsStatsVals),3);
            for k = 1:length(IndsStatsVals)
                XYZmmTmp{1}(k,:) = XYZvoxTmp{1}(k,:)*Vol.mat(1:3,1:3)+Vol.mat(1:3,4)';
            end
        end
    end
end

%negative stats vals
if(any(Data3D(:)<Thresh(2)))
    StatsValsTmp{2} = Data3D(Data3D(:)<Thresh(2)); %statistics values NEGATIVE
    IndsStatsVals   = find(Data3D(:)<Thresh(2));   %corresponding indices
    XYZvoxTmp{2}    = zeros(length(IndsStatsVals),3); %nVoxels-x-3 (spatial dimensions)
    [XYZvoxTmp{2}(:,1),XYZvoxTmp{2}(:,2),XYZvoxTmp{2}(:,3)] = ind2sub(Vol.dim,IndsStatsVals); %indices to 3D/subscript
    
    if(nargout>2) %Calculate only if mm-Coords are needed.
        %do trafo from voxels to world (mm-coords)
        TestDiff = Vol.mat(1:3,1:3)-diag(diag(Vol.mat(1:3,1:3))); %if this contains only zeros then we can do this very fast!
        if(all(TestDiff(:)==0)) %rotation&scaling submatrix is diagonal --> we can do it fast with bsxfun
            XYZmmTmp{2}=bsxfun(@plus,bsxfun(@times,XYZvoxTmp{2},diag(Vol.mat(1:3,1:3))'),Vol.mat(1:3,4)'); %fast affine transformation using bsxfun (3D submatrix is diagonal)
        else %not diagonal
            disp('WARNING(neg): VOXEL2WORLD 3D submatrix is NOT diagonal!!! (NB: usually any normalized (MNI) space should be diagonal in the 3D subspace.)'); %usually any normalized (MNI) space should be diagonal in the 3D subspace.
            XYZmmTmp{2}=zeros(length(IndsStatsVals),3);
            for k = 1:length(IndsStatsVals)
                XYZmmTmp{2}(k,:) = XYZvoxTmp{2}(k,:)*Vol.mat(1:3,1:3)+Vol.mat(1:3,4)';
            end
        end
    end
end

%% reform inputs
StatsVals = [StatsValsTmp{1}; StatsValsTmp{2}];
XYZvox    = [XYZvoxTmp{1};    XYZvoxTmp{2}];
XYZmm     = [XYZmmTmp{1};     XYZmmTmp{2}];

end