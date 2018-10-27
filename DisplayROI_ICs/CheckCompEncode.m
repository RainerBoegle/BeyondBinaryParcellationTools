function [H,ax] = CheckCompEncode(CompEncodeDir)
% This function will "look" through the results of the component encode process 
% and check for components that are in the UniqueMap and those that are in the OverlapMap
% such that one can see which components these are, and if there are any that are purely Unique
% or purely Overlap, ie not in the UniqueMap.
%
% This will just be a basic plot nothing too sophisticated. ;)
%
%Usage:
%       H = CheckCompEncode(CompEncodeDir);
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (28.08.2015): initial implementation

%% check input & load files if possible
%CompEncodeDir --> get nifti-files
try
    if(~exist(CompEncodeDir))
        CompEncodeDir = spm_select(1,'dir','Select directory with component encode NIFTI-files...');
    else
        %search for files
        [Vols,Files,NIFTI_files] = GetCompEncodeFiles(CompEncodeDir);
    end
catch CATCH_CompEncodeDir
    disp_catch(CATCH_CompEncodeDir,[mfilename,'>CATCH_CompEncodeDir'],'CATCH_CompEncodeDir'); %assign & display error as warning
    CompEncodeDir = spm_select(1,'dir','Select directory with component encode NIFTI-files...');
    %search for files
    [Vols,Files,NIFTI_files] = GetCompEncodeFiles(CompEncodeDir);
end

%% get number of components
Y_Base2Encode = Vols{find(~cellfun(@isempty,strfind(NIFTI_files,[filesep,'Base2Encode'])))}.private.dat(:);
NComp = length(dec2bin(max(Y_Base2Encode(:)))); %the maximum will either be the component number in binary and because of Base2Encode the length of the binary number (which is repressented as string) or maximally the mixture of all components; -therefore the length still indicates the number of components.

%% examine contents of UniqueMap 
Y_UniqueOnly = Vols{find(~cellfun(@isempty,strfind(NIFTI_files,[filesep,'UniqueOnly'])))}.private.dat(:);  

CompIn_UniqueOnlyVol = zeros(NComp,1);
if(~isempty(find(Y_UniqueOnly~=0)))
    Inds = unique(Y_UniqueOnly(Y_UniqueOnly~=0));
    if(any(Inds<=0))
        disp('WARNING there are zero and/or negative indices in UniqueOnlyMap??? Something went wrong, please check. Will remove those for plotting.');
        Inds(Inds<=0) = [];
    end
    CompIn_UniqueOnlyVol(Inds) = 1;
end

%% examine contents of OverlapOnly Volume via contents of WinnerTakeAllEncode (i.e. Mask with OverlapMap and see which are there...)
Y_OverlapOnly         = Vols{find(~cellfun(@isempty,strfind(NIFTI_files,[filesep,'OverlapOnly'])))}.private.dat(:);  
Y_WinnerTakeAllEncode = Vols{find(~cellfun(@isempty,strfind(NIFTI_files,[filesep,'WinnerTakeAllEncode'])))}.private.dat(:);  

CompIn_OverlapOnlyVol_fromWinnerTakeAllEncode = zeros(NComp,1);
IndsUnique = unique(Y_WinnerTakeAllEncode(Y_OverlapOnly~=0));
IndsUnique(IndsUnique==0) = [];
CompIn_OverlapOnlyVol_fromWinnerTakeAllEncode(IndsUnique) = 1;

%% and using Base2Encode
UniqueBase2CompMixIn_OverlapOnlyVol_fromBase2Encode = unique(Y_Base2Encode(Y_OverlapOnly~=0));
CompIn_OverlapOnlyVol_fromBase2Encode = zeros(NComp,1);
for Ind = 1:length(UniqueBase2CompMixIn_OverlapOnlyVol_fromBase2Encode)
    tempBinary = arrayfun(@str2num,dec2bin(UniqueBase2CompMixIn_OverlapOnlyVol_fromBase2Encode(Ind))); %take current component mix number (base2encode) make binary and convert to a vector of 0&1
    tempBinary = tempBinary(end:-1:1); %reverse such that order is right
    CompIn_OverlapOnlyVol_fromBase2Encode(tempBinary~=0) = 1; %assign those positions with a 1.
end
    

%% display results
H{1} = figure;
ax{1}(1) = subplot(3,1,1); bar((1:NComp)',CompIn_UniqueOnlyVol); title('Components contained in UniqueOnly Volume...');
ax{1}(2) = subplot(3,1,2); bar((1:NComp)',CompIn_OverlapOnlyVol_fromWinnerTakeAllEncode); title('Components contained in OverlapOnly Volume [based on WinnerTakeAllEncode Volume]...');
ax{1}(3) = subplot(3,1,3); bar((1:NComp)',CompIn_OverlapOnlyVol_fromBase2Encode); title('Components contained in OverlapOnly Volume [based on Base2Encode Volume]...');

%% display trees of Base2Encode
IndsBase2Enc = unique(Y_Base2Encode(Y_Base2Encode~=0));
if(any(IndsBase2Enc<=0))
    disp('WARNING there are zero and/or negative indices in Y_Base2Encode??? Something went wrong, please check. Will remove those for plotting.');
    IndsBase2Enc(IndsBase2Enc<=0) = [];
end
BinaryMatrix_UniqueInds_Base2Encode = arrayfun(@str2num,dec2bin(IndsBase2Enc));
BinaryMatrix_UniqueInds_Base2Encode = BinaryMatrix_UniqueInds_Base2Encode(:,end:-1:1); %turn around to move first bit, i.e. IC1 to first index
BinaryMatrix_UniqueInds_Base2Encode = BinaryMatrix_UniqueInds_Base2Encode'; %flip rows & columns
Overlaps_BinaryMatrix_UniqueInds_Base2Encode = BinaryMatrix_UniqueInds_Base2Encode;
DelInds = []; %init indices of pure connections that are removed because "they are not a tree", i.e. have no connections
for Ind = 1:size(Overlaps_BinaryMatrix_UniqueInds_Base2Encode,2)
    if(sum(Overlaps_BinaryMatrix_UniqueInds_Base2Encode(:,Ind))==1)
        DelInds = [DelInds Ind];
    end
end
Overlaps_BinaryMatrix_UniqueInds_Base2Encode(:,DelInds)=[];
H{2} = figure; 
subplot(1,2,1); imagesc(Overlaps_BinaryMatrix_UniqueInds_Base2Encode); title('IC connections for all possible overlaps that were found in Base2Encode (i.e. all indices that are not a power of 2 converted to binary).'); ylabel('ICs'); axis('xy') %usually better to view it this way.
subplot(1,2,2); imagesc(sum(Overlaps_BinaryMatrix_UniqueInds_Base2Encode,2)); title('IC connections for all possible overlaps that were found in Base2Encode (i.e. all indices that are not a power of 2 converted to binary).'); ylabel('ICs'); axis('xy') %usually better to view it this way.
colormap('hot'); colorbar;

H{3} = figure; dendrogram(linkage(Overlaps_BinaryMatrix_UniqueInds_Base2Encode)); title('simple dendrogram of clusters');

%% done.
disp('Done checking Component Encode results.');


end


%% subfunction for laoding data
function [Vols,Files,NIFTI_files] = GetCompEncodeFiles(CompEncodeDir)
% get NIFTI-files from CompEncodeDir

%try to find these NIFTI-files
% Base2Encode.nii     %Not for display
% Log2Base2Encode.nii %Display special colors
% SumMap.nii          %Display with counter colorbar
% WinnerTakeAllEncode %Display
% SimpleEncode.nii    %Display
% OverlapOnly.nii     %Display
% UniqueOnly.nii      %Display
Files = cell(7,1);
Files{1} = 'Base2Encode.nii';
Files{2} = 'Log2Base2Encode.nii';
Files{3} = 'SumMap.nii';
Files{4} = 'WinnerTakeAllEncode.nii';
Files{5} = 'SimpleEncode.nii';
Files{6} = 'OverlapOnly.nii';
Files{7} = 'UniqueOnly.nii';

NIFTI_files    = cell(7,1);
NIFTI_files{1} = [CompEncodeDir,filesep,'Base2Encode.nii'];
NIFTI_files{2} = [CompEncodeDir,filesep,'Log2Base2Encode.nii'];
NIFTI_files{3} = [CompEncodeDir,filesep,'SumMap.nii'];
NIFTI_files{4} = [CompEncodeDir,filesep,'WinnerTakeAllEncode.nii'];
NIFTI_files{5} = [CompEncodeDir,filesep,'SimpleEncode.nii'];
NIFTI_files{6} = [CompEncodeDir,filesep,'OverlapOnly.nii'];
NIFTI_files{7} = [CompEncodeDir,filesep,'UniqueOnly.nii'];

Vols = cell(length(Files),1);
for IndFile = 1:length(Files)
    if(~exist(NIFTI_files{IndFile},'file'))
        error(['Can not find "',Files{IndFile},'" in "',CompEncodeDir,'".']);
    else
        Vols{IndFile} = spm_vol(NIFTI_files{IndFile});
    end
end

end