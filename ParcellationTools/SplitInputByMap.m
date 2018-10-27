function [InfoStruct] = SplitInputByMap(InputMapPath,SplitMapPath,Limits,OutputDir)
% This function allows to split an Input image in multiple images that is
% masked by the possible overlaps in the SumMap.
% I.e. the sum map is turned into a mask from it's lowest overlap value 1 == uniques
% to its highest value.
%
%Inputs:
%        InputMapPath    (string/char)      Indicates the path to the NIFTI-file that should be split up according to the Map.
%        SplitMapPath    (string/char)      Indicates the path to the NIFTI-file that repressents the Map that is then used for mask creation together with "Limits" input for splitting the Input 
%        Limits          (N-x-2 matrix)     Each row indicates the lower (1st Column) and upper value (2nd Column) that should be used to make a mask out of SplitMapPath (i.e. these values need to be appropriate for SplitMap).
%                        (or string/char)    EXAMPLE: 1.For the Input [1, 1; 2, 2; 3, 3; 4, inf]
%                                                      This mean there will be a split for the values in SplitMap
%                                                      that are == 1 and those == 2 and those == 3, as well as for all >= 4.
%                                            EXAMPLE: 2.For the string input 'Unique' instead of the matrix input to get a split for each unique value in SplitMap (except ZERO). 
%                                                       E.g. if SplitMap contains values from 1 to 6, this will be the same as supplying the input matrix [1, 1; 2, 2; 3, 3; 4, 4; 5, 5; 6, 6]. 
%                                            WARNING: if SplitMap is not discrete this will produce strange results!
%        OutputDir       (string/char)      Indicates the path to the which the resulting split NIFTI-files should be written.                           
%                     [empty/not given]     If empty or not given, then a directory will be created in the directory of InputMapPath for outputting the splits.
%
%
%Usage:
%       [InfoStruct] = SplitInputByMap(InputMapPath,SplitMapPath,Limits,OutputDir);
%       [InfoStruct] = SplitInputByMap(InputMapPath,SplitMapPath,'Unique');
%       [InfoStruct] = SplitInputByMap(InputMapPath,SplitMapPath,[1, 1; 2, 2; 3, 3; 4, inf]);
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (28.10.2016): initial implementation based on previous tests

%% Check inputs
if(~exist('InputMapPath','var'))
    InputMapPath = spm_select(1,'image','Select InputMap that should be split up...');
else
    if(isempty(InputMapPath))
        InputMapPath = spm_select(1,'image','Select InputMap that should be split up...');
    else
        if(~ischar(InputMapPath))
            if(iscellstr(InputMapPath))
                InputMapPath = InputMapPath{1}; %taking first input
            else
                error('InputMapPath must be a string/char-array pointing to the Map that should be split up.');
            end
        end
    end
end
disp(['Will split up InputMap "',InputMapPath,'"...']);

if(~exist('SplitMapPath','var'))
    SplitMapPath = spm_select(1,'image','Select SplitMap that is the template for creating the split masks based on Limits...');
else
    if(isempty(SplitMapPath))
        SplitMapPath = spm_select(1,'image','Select SplitMap that is the template for creating the split masks based on Limits...');
    else
        if(~ischar(SplitMapPath))
            if(iscellstr(SplitMapPath))
                SplitMapPath = SplitMapPath{1}; %taking first input
            else
                error('SplitMapPath must be a string/char-array pointing to the Map that provides the split mask based on Limits that are applied.');
            end
        end
    end
end

if(~exist('Limits','var'))
    error('Limits not input, can be a Nsplits-x-2 matrix indicating lower and upper limits for the SplitMap OR the string ''Unique'' which will split each unique value of SplitMap.');
else
    if(isempty(Limits))
        error('Limits is empty! Limits can be a Nsplits-x-2 matrix indicating lower and upper limits for the SplitMap OR the string ''Unique'' which will split each unique value of SplitMap.');
    else
        if(isnumeric(Limits))
            if(size(Limits,2)~=2)
                error('Limits must be a Nsplits-x-2 matrix indicating lower and upper limits for the SplitMap OR the string ''Unique'' which will split each unique value of SplitMap.');
            end
        else
            if(ischar(Limits))
                if(~strcmpi(Limits,'Unique'))
                    error('Limits must be the string ''Unique'' which will split each unique value of SplitMap OR a Nsplits-x-2 matrix indicating lower and upper limits for the SplitMap.');
                end
            else
                error('Limits must be the string ''Unique'' which will split each unique value of SplitMap OR a Nsplits-x-2 matrix indicating lower and upper limits for the SplitMap.');
            end
        end
    end
end

if(~exist('OutputDir','var'))
    OutputDir = []; %fill later
else
    if(~isempty(OutputDir))
        if(~exist(OutputDir,'dir'))
            disp(['OutputDir "',OutputDir,'" does not exist, will create it...']);
            mkdir(OutputDir);
        else
            disp(['Will save splits into directory "',OutputDir,'"...']);
        end
    end
end

%% Setup OutputNameBase
[~,FNameInput] = fileparts(InputMapPath);
[~,FNameMap]   = fileparts(SplitMapPath); %this allows other applications besides SumMap.
OutputNameBase = [FNameInput,'_SplitBy_',FNameMap];

%% check output dir
if(isempty(OutputDir)) %create in basedir of InputMapPath
    BaseDir = fileparts(InputMapPath);
    OutputDir = [BaseDir,filesep,OutputNameBase];
    disp(['OutputDir "',OutputDir,'" does not exist, will create it...']);
    mkdir(OutputDir);
end

%% Check Limits
Limits = CheckLimits(SplitMapPath,Limits);

%% Iterate over Limits and create matlabbatch for using ImCalc to create split files
matlabbatch = cell(size(Limits,1),1);
for Ind = 1:size(Limits,1)
    LowLimit = Limits(Ind,1);
    UpLimit  = Limits(Ind,2);
    
    %% setup matlabbatch
    OutputName       = [OutputNameBase,'_',num2str(LowLimit),num2str(UpLimit),'.nii'];
    CurrMBatch       = SetupBatch(InputMapPath,SplitMapPath,LowLimit,UpLimit,OutputDir,OutputName);
    matlabbatch{Ind} = CurrMBatch{1};
end

%% apply ImCalc
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

%% save InfoStruct
InfoStruct = struct('InputMapPath',InputMapPath,'SplitMapPath',SplitMapPath,'Limits',Limits,'OutputDir',OutputDir,'matlabbatch',{matlabbatch});
disp(['Saving InfoStruct to "',OutputDir,filesep,'InfoStruct.mat".']);
save([OutputDir,filesep,'InfoStruct.mat'],'InfoStruct');

%% done.
disp(' ');
disp('Done.');
disp(' ');


end

%% subfunctions

function Limits = CheckLimits(SplitMapPath,Limits)
% setup Limits

if(isnumeric(Limits))
    return; %this should always work i.e. at least execute, but of course it depends on the user to specify it correctly, such that the result is correct.
else
    if(strcmpi(Limits,'Unique'))
        %get from volume
        Vol = spm_vol(SplitMapPath);
        Data= Vol.private.dat(:,:,:,Vol.n(1));
        Uniques = unique(Data(Data(:)~=0));
        %check if integer
        if(any(mod(Uniques,1)~=0))
            disp('WARNING: SplitMap "',SplitMapPath,'" does not only contain integers/whole numbers!!! Using the ''Unique'' option might lead to strange results!');
            disp('WARNING: SplitMap "',SplitMapPath,'" does not only contain integers/whole numbers!!! Using the ''Unique'' option might lead to strange results!');
            disp('WARNING: SplitMap "',SplitMapPath,'" does not only contain integers/whole numbers!!! Using the ''Unique'' option might lead to strange results!');
        end
        Limits = zeros(length(Uniques),2);
        for Ind = 1:length(Uniques)
            Limits(Ind,:) = [Uniques(Ind), Uniques(Ind)];
        end
    end
end

end



function [matlabbatch] = SetupBatch(InputMapPath,SplitMapPath,LowLimit,UpLimit,OutputDir,OutputName)
%Setup batch

matlabbatch{1}.spm.util.imcalc.input          = cell(2,1);
matlabbatch{1}.spm.util.imcalc.input{1}       = InputMapPath;
matlabbatch{1}.spm.util.imcalc.input{2}       = SplitMapPath;
matlabbatch{1}.spm.util.imcalc.output         = OutputName; %'WinnerTakeAll_SumMapMask_1.nii';
matlabbatch{1}.spm.util.imcalc.outdir{1}      = OutputDir;
matlabbatch{1}.spm.util.imcalc.expression     = 'i1.*((i2>=LowLimit).*(i2<=UpLimit))';
matlabbatch{1}.spm.util.imcalc.var(1).name    = 'LowLimit';
matlabbatch{1}.spm.util.imcalc.var(1).value   = LowLimit;
matlabbatch{1}.spm.util.imcalc.var(2).name    = 'UpLimit';
matlabbatch{1}.spm.util.imcalc.var(2).value   = UpLimit;
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype  = 16;

end