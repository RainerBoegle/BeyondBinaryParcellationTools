function  [InfoStruct,NIIOutputPaths,InfoStructPath] = genICBase2EncThreshOverlap(ResultsStruct,pThresh,OutputDir)
% This function uses the ResultsStruct (from checkICsWithMasks.m) and a percentage Threshold [0 1) to generate 
% Base2Encode Overlap Maps & SumMaps PER IC for all other ICs that overlap it more than pThresh,
% as determined from ResultsStruct MaxScaledSumContributions and DiceCoefficients.
%
% Three different Versions are saved INCLUDING a description in a InfoStruct that is saved in a *.mat-file.
% All maps are saved as 4D NIFTIs.
% For the Overlap of both and the sqrt(product of both types) to create the final score that is then thresholded. 
% NB: this will reduce the number of overlapping components.
%
%
% THE MAIN IDEA behind this overlap plotting is similar to the Base2Encode Map,
% but here it is contained for a single component AT A TIME
% and due to the thresholding it is for a reduced number of overlaps
% that start counting for each IC i.e. the numbers don't get too high.
%
% Steps for each kind of Map after thresholding:
% 1.for each IC make this IC into a Mask
% 2.add each overlapping IC as powers of 2
% 3.save this in 4D array
% when done save this as 4D NIFTI with name and threshold. 
%
%
%Usage:
%       [InfoStruct,NIIOutputPaths,InfoStructPath] = genICBase2EncThreshOverlap(ResultsStruct,pThresh,OutputDir);
%       [InfoStruct,NIIOutputPaths,InfoStructPath] = genICBase2EncThreshOverlap(ResultsStruct,0.1);   %Threshold at 10%(default) contribution and select OutputDir manually
%       [InfoStruct,NIIOutputPaths,InfoStructPath] = genICBase2EncThreshOverlap(ResultsStruct);       %SAME AS ABOVE Threshold at 10%(default) contribution and select OutputDir manually
%       [InfoStruct,NIIOutputPaths,InfoStructPath] = genICBase2EncThreshOverlap(checkICsWithMasks()); %generate a ResultsStruct & then SAME AS ABOVE. Threshold at 10%(default) contribution and select OutputDir manually
%       [InfoStruct,NIIOutputPaths,InfoStructPath] = genICBase2EncThreshOverlap(checkICsWithMasks(),pThresh,OutputDir); 
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0(11.10.2016): initial implementation based on previous tests


%% check inputs
if(~isstruct(ResultsStruct))
    error('Input ResultsStruct is not a struct!');
end

if(~exist('pThresh','var'))
    pThresh = 0.1;
else
    if(isempty(pThresh))
        pThresh = 0.1;
    else
        if(pThresh>=1||pThresh<0)
            error('pThresh must be between 0 (inclusive) and 1 (exclusive), i.e. [0 1) or 0<= pThresh <1.');
        end
    end
end
disp(['Using pThresh= ',num2str(pThresh),'...']);

if(~exist('OutputDir','var'))
    OutputDir = spm_select(1,'dir','Select Output Directory for Results...');
else
    if(isempty(OutputDir))
        OutputDir = spm_select(1,'dir','Select Output Directory for Results...');
    else
        if(~ischar(OutputDir))
            if(iscellstr(OutputDir)&&length(OutputDir)==1)
                OutputDir = OutputDir{1};
            else
                error('OutputDir must be a string/cellstring that indicates the directory that should be used for storing the Results.');
            end
        end
    end
end
if(~exist(OutputDir,'dir'))
    disp(['Directory "',OutputDir,'" does not exist, will create it...']);
    mkdir(OutputDir);
end


%% get data & Threshold
InfoStruct.ThreshNIIpaths = ResultsStruct.ThreshNIIpaths;
InfoStruct.NICs           = length(InfoStruct.ThreshNIIpaths);
InfoStruct.VICs           = spm_vol(InfoStruct.ThreshNIIpaths);
InfoStruct.pThresh        = pThresh;
InfoStruct.Data           = cell(3,1);
InfoStruct.DataName       = cell(3,1);
InfoStruct.ThreshData     = cell(3,1);

%MaxScaledSumContributions
InfoStruct.DataName{1} = 'SumContributions';
InfoStruct.DataMat{1}  = ResultsStruct.MaxScaleSumContributions; %using max scaled one for simplicity

%DiceCoeff
InfoStruct.DataName{2} = 'DiceCoeff';
InfoStruct.DataMat{2}  = ResultsStruct.DiceCoeff;

%ProductCoeff
InfoStruct.DataName{3} = 'ProductCoeff';
InfoStruct.DataMat{3}  = ResultsStruct.ProductCoeff;

%% for each data create output
NIIOutputPaths = cell(length(InfoStruct.DataName),2);
for IndData = 1:length(InfoStruct.DataName)
    disp(' ');
    disp(['Creating individual Base2Encodes per IC using data from "',InfoStruct.DataName{IndData},'"...']);
    %% create 4D Data in Array
    Vtmp   = InfoStruct.VICs{1}; %use this as a template
    Base2Enc4D = zeros(Vtmp.dim(1),Vtmp.dim(2),Vtmp.dim(3),InfoStruct.NICs);
    SumMap4D   = zeros(Vtmp.dim(1),Vtmp.dim(2),Vtmp.dim(3),InfoStruct.NICs);
    for IndIC = 1:InfoStruct.NICs
        Mask3D = double(InfoStruct.VICs{IndIC}.private.dat(:,:,:,InfoStruct.VICs{IndIC}.n(1))~=0);
        [Base2Enc4D(:,:,:,IndIC),SumMap4D(:,:,:,IndIC),OK] = createBase2Enc_n_SumMap(InfoStruct.VICs,Mask3D,InfoStruct.DataMat{1}(IndIC,:),pThresh);
        if(~OK)
            disp(['The problem occurred at IC ',num2str(IndIC),'... check later...']);
        end
    end
    disp('...done');
    
    %% write out Results to OutputDir as NIFTI
    disp(['Writing out Results to Directory "',OutputDir,'" as 4D-NIFTI...']);
    %Base2Enc4D
    NIIOutputPaths{IndData,1} = [OutputDir,filesep,'IndividualBase2Enc_',InfoStruct.DataName{IndData},'_pThresh',num2str(pThresh*100),'.nii'];
    Write4DNIFTI(Vtmp,Base2Enc4D,NIIOutputPaths{IndData,1});
    
    %SumMap4D
    NIIOutputPaths{IndData,2} = [OutputDir,filesep,'IndividualSumMap_',InfoStruct.DataName{IndData},'_pThresh',num2str(pThresh*100),'.nii'];
    Write4DNIFTI(Vtmp,SumMap4D,NIIOutputPaths{IndData,2});
    disp('...done.');
end
disp(' ');

%% save InfoStruct to OutputDir as well.
disp(['Saving InfoStruct to Directory "',OutputDir,'"...']);
InfoStructPath = [OutputDir,filesep,'IndividualBase2Enc_InfoStruct.mat'];
save(InfoStructPath,'InfoStruct');
disp('...done.');

%% Done.
disp(' ');
disp('DONE.');
disp(' ');

end

%% subfunctions
function [Base2Enc,SumMap,OK] = createBase2Enc_n_SumMap(Vols,Mask3D,InputData,pThresh)
% This function does the Base2Encode & uses the inital IC as a Mask and the InputData is the overlap
% with the other ICs which is sorted and then thresholded using pThresh 
% AND only those volumes are used that survived the theshold and these are added onto the mask of
% that ICs (given by Mask3D) 
%
% --> Then we get the overlaps in that specific IC as Base2Encode. Done.

Dims = size(Mask3D);

%% Init Base2Enc & SumMap
Base2Enc = zeros(Dims(1),Dims(2),Dims(3));
SumMap   = zeros(Dims(1),Dims(2),Dims(3));
OK      = 0; %init not okay. (empty/zero data output)

%% determine which Inputs should be collected
[SortedVals,SortInds] = sort(InputData,'descend');
SortInds = SortInds(SortedVals>pThresh); %pick those that are appropriate.
if(isempty(SortInds))
    disp('SortInds is empty, i.e. no overlaps with this Input???');
    return;
end

%% do Base2Encode
for Ind = 1:length(SortInds)
    Base2Enc = Base2Enc + (2^(Ind-1)).*double((Vols{SortInds(Ind)}.private.dat(:,:,:,Vols{SortInds(Ind)}.n(1))~=0).*Mask3D);
    SumMap   = SumMap   +              double((Vols{SortInds(Ind)}.private.dat(:,:,:,Vols{SortInds(Ind)}.n(1))~=0).*Mask3D);
end

%% at this point we should be fine
OK = 1;

end

%% write NIFTI
function [Vout] = Write4DNIFTI(Vtmp,Data4D,OutputPath)
% Write out the 4D NIFTI

%% settings
Vout = rmfield(Vtmp,'private');
Vout.fname = OutputPath;
if(Vout.dt(1)<16)
    Vout.dt(1) = 16; %just to be sure make encoding depth higher than needed.
end

%% write data
for IndIC = 1:size(Data4D,4)
    Vout.n(1) = IndIC;
    spm_write_vol(Vout,Data4D(:,:,:,IndIC));
end

%% get structures of output
Vout = spm_vol(OutputPath); %reload is the best way of getting all the structures

end