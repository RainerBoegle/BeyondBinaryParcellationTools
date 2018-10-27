function [OutputPath,InfoStruct] = CreateSpecialEncodeMap(CompEncodeDir,CombineOverlapNumbers)
% This function will create a special encoding map that combines the
% properties of the SumMap and the OverlapsOnlyMap and the UniqueOnlyMap.
%
% Basic Idea:
% The user can choose the number of overlaps to be assigned to a certain color
% and all the uniques get assigned to colors different from these colors.
% This way the user can color 2-way overlaps differently from 3-way to N-way overlaps. 
% E.g.: dark-blue 2-way overlap and 3-way overlaps in blue and overything higher in light blue or cyan. 
% Then add all ICs in different colors on top. 
% Then the user get's an information figure regarding the colors and then the overlay.
%
% Hopefully this will help with the display of uniques and overlays,
% because we need to know a little bit more than just where are the uniques
% and where the overlays, but how many overlaps occur.
% e.g. 
%     CombineOverlapNumbers = {[2,3];[4,inf]};     %2&3 are one overlap class and 4:maximum overlap are one class.
%     CombineOverlapNumbers = {[2,3];[4];[5,inf]}; %2&3 are one overlap class and 4 is another and 5:maximum overlap are one class.
%
%
%Usage:
%       [OutputPath,InfoStruct] = CreateSpecialEncodeMap(CompEncodeDir,CombineOverlapNumbers);
%
%
%V1.1
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.1: (26.10.2016): added possibility to define CombineOverlapNumbers as {[2,3];[4,inf]} as indicating that all numbers above 3 are the last combined overlap. V1.0: (06.10.2016): initial implementation based on previous tests

%% Check Inputs
if(~exist('CompEncodeDir','var'))
    CompEncodeDir = spm_select(1,'dir','Select Component Encode Directory...');
else
    if(isempty(CompEncodeDir))
        CompEncodeDir = spm_select(1,'dir','Select Component Encode Directory...');
    elseif(~ischar(CompEncodeDir))
        if(iscellstr(CompEncodeDir))
            CompEncodeDir = CompEncodeDir{1}; %taking the first one
        else
            error('CompEncodeDir must be string indicating the path to the Component Encode Directory...');
        end
    end
end
if(~exist(CompEncodeDir,'dir'))
    error(['Folder "',CompEncodeDir,'" could not be found!']);
end

if(~exist('CombineOverlapNumbers','var'))
    error('CombineOverlapNumbers not given!');
else
    if(~iscell(CombineOverlapNumbers))
        error('ColorCombineOverlapNumbers must be a cell!');
    end    
end

%% get files from Component Encode Directory & get  DataStruct, the volume structures(SPM) for SumMap and UniqueOnly
disp('Loading data...');
DataStructPath = spm_select('FPList',CompEncodeDir,'^DataStruct.mat');
if(isempty(DataStructPath))
    error(['Could not find DataStruct.mat in directory "',CompEncodeDir,'"!']);
else
    InfoStruct = load(DataStructPath);
end

SumMapPath = spm_select('FPList',CompEncodeDir,'^SumMap.nii');
if(isempty(SumMapPath))
    error(['Could not find SumMap.nii in directory "',CompEncodeDir,'"!']);
else
    InfoStruct.V_SumMap = spm_vol(SumMapPath);
end

UniqueOnlyPath = spm_select('FPList',CompEncodeDir,'^UniqueOnly.nii');
if(isempty(UniqueOnlyPath))
    error(['Could not find UniqueOnly.nii in directory "',CompEncodeDir,'"!']);
else
    InfoStruct.V_UniqueOnly = spm_vol(UniqueOnlyPath);
end

%% get data in 3D
YSumMap     = InfoStruct.V_SumMap.private.dat(:,:,:,InfoStruct.V_SumMap.n(1));
YUniqueOnly = InfoStruct.V_UniqueOnly.private.dat(:,:,:,InfoStruct.V_UniqueOnly.n(1));

%% check CombineOverlapNumbers in case inf was used 
TestInf = cellfun(@sum,cellfun(@isinf,CombineOverlapNumbers,'UniformOutput',false));
if(any(TestInf~=0))
    if(length(find(TestInf~=0))>1)
        error('There can only be one "inf" case! (and it has to be the last entry)');
    else
        if(TestInf(end)~=1)
            error('There can only be one "inf" case and it has to be the last entry!');
        end
    end
    FinalOverlap = max(YSumMap(:));
    if(CombineOverlapNumbers{end}(~isinf(CombineOverlapNumbers{end}))<=FinalOverlap)
        CombineOverlapNumbers{end} = CombineOverlapNumbers{end}(~isinf(CombineOverlapNumbers{end})):FinalOverlap;
    else
        %remove this overlap category
        CombineOverlapNumbers(end) = [];
    end
end

%% Create Special Encode Step 1: encode overlaps by using settings in CombineOverlapNumbers
disp('Create Special Encode Step 1: encode overlaps by using settings in CombineOverlapNumbers...');
YSpecialEncode = zeros(size(YSumMap));
for IndAssign = 1:length(CombineOverlapNumbers)
    CurrCombinedOverlaps = CombineOverlapNumbers{IndAssign};
    for IndCurr = 1:length(CurrCombinedOverlaps)
        YSpecialEncode(YSumMap==CurrCombinedOverlaps(IndCurr)) = IndAssign;
    end
end

%% Create Special Encode Step 2: add UniqueOnly values "on top"
disp('Create Special Encode Step 2: add UniqueOnly values "on top"...');
ICsUnique = unique(YUniqueOnly);
ICsUnique(ICsUnique==0) = []; %remove zero
for IndICs = 1:length(ICsUnique)
    YSpecialEncode(YUniqueOnly==ICsUnique(IndICs)) = length(CombineOverlapNumbers) + IndICs;
end
disp('...done.');

%% save NIFTI
disp(['Saving SpecialEncodeMap.nii to directory "',CompEncodeDir,'".']);
OutputPath = [CompEncodeDir,filesep,'SpecialEncodeMap.nii'];
Vout       = rmfield(InfoStruct.V_UniqueOnly,'private');
Vout.fname = OutputPath;
InfoStruct.V_SpecialEncode = spm_write_vol(Vout,YSpecialEncode);

%% set InfoStruct
InfoStruct.CombineOverlapNumbers = CombineOverlapNumbers;
InfoStruct.ICcorrespondence      = [(length(CombineOverlapNumbers)+(1:length(ICsUnique))'),ICsUnique(:)];
InfoStruct.ICcorrespondenceInfo  = 'First column is the numbers in the SpecialEncodeMap and the second column is the corresponding number of the original ICs.';

%% save InfoStruct
disp(['Saving InfoStruct to SpecialEncodeInfoStruct.mat in directory "',CompEncodeDir,'".']);
save([CompEncodeDir,filesep,'SpecialEncodeInfoStruct.mat'],'InfoStruct');

%% Done.
disp(' ');
disp('DONE.');
disp(' ');

end
