function [OutputDir,CompEncodeOutPaths,OverlapMatStruct,DataStruct] = CreateCompEncodeFiles(InputFiles,ThreshType,Thresh,OutputDir)
% This function can be used to create the Component Encode Maps/NIFTI-files and Overlap Matrix,
% based on the input files and an optional threshold.
% The results will then be output to the specified directory "OutputDir".
%
%Usage:
%       [OutputDir,OverlapMatStruct] = CreateCompEncodeFiles(InputFiles,ThreshType,Thresh,OutputDir);
%       [OutputDir,OverlapMatStruct] = CreateCompEncodeFiles(InputFiles,ThreshType,Thresh);           %select output directory       
%       [OutputDir,OverlapMatStruct] = CreateCompEncodeFiles(InputFiles,ThreshType,    [],OutputDir); %do not apply a threshold, i.e. all non-zero values are used.
%       [OutputDir,OverlapMatStruct] = CreateCompEncodeFiles(InputFiles,       'P',   0.5,OutputDir); %take top 50% of voxels 
%       [OutputDir,OverlapMatStruct] = CreateCompEncodeFiles(InputFiles,   'Stats',[2 -4],OutputDir); %take voxels above stats value 2(a.u.) and those below stats value -4(a.u.).
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (22.03.2016): initial implementation based on earlier test versions 

%% check inputs
%InputFiles?
try
    if(~iscellstr(InputFiles))
        if(ischar(InputFiles))
            InputFiles = cellstr(InputFiles);
        end
    end
catch CATCH_InputFiles
    disp_catch(CATCH_InputFiles,'CreateOverlapMat>Input:InputFiles','CATCH_InputFiles');
    InputFiles = cellstr(spm_select(Inf,'image','Select thresholded maps for creation of OverlapMat...'));
end

%ICnums?
[ok,ICnums] = findICnums(InputFiles); %get IC nums assuming they are listed at the end of the file... (subfunction below)
if(ok) %could extract ICnums
    if(~issorted(ICnums)) %are they sorted?
        disp('IC# were found, BUT NOT SORTED! Will sort everything in ascending order...');
        [ICnums,sortInds] = sort(ICnums);
        InputFiles = InputFiles(sortInds); %InputFiles must be a cellstr otherwise this doesn't work therefore the hassle above when checking the inputs.
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
        case 'Stats'
            ThreshType  = 'Stats'; %just to be save.
            ThreshStats = Thresh;
        case 'HalfMax'
            ThreshType  = 'HalfMax'; %just to be save.
            ThreshStats = Thresh;
        otherwise
            error(['Unknown ThreshType: "',ThreshType,'"?!']);
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
        disp([ThreshType,': Using threshold type statistics values (a.u.) above and below threshold.']);          
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
   
%OutputDir
if(~exist('OutputDir','var'))
    OutputDir = uigetdir(pwd,'Select directory to save results of component encoding...');
else
    if(isempty(OutputDir))
        OutputDir = uigetdir(pwd,'Select directory to save results of component encoding...');
    elseif(ischar(OutputDir))
        if(isdir(OutputDir))
            disp(['Will save results to directory "',OutputDir,'"...']);
        else
            mkdir(OutputDir);
        end
    elseif(iscellstr(OutputDir))
        OutputDirTmp = OutputDir; clear OutputDir
        OutputDir    = OutputDirTmp{1}; clear OutputDirTmp
        if(isdir(OutputDir))
            disp(['Will save results to directory "',OutputDir,'"...']);
        else
            mkdir(OutputDir);
        end
    else
        error('"OutputDir" must a char/string or cellstring!');
    end
end

%% do calculations
%% extract thresholds and create SimpleEncode, Base2Encode, OverlapOnly and UniqueOnly volume
CompEncodeOutPaths= cell(7,1); %SimpleEncode; Base2Encode; OverlapOnly; UniqueEncode

%init using first selected IC
%SimpleEncode
V_SimpleEnc = spm_vol(InputFiles{1});
Y_SimpleEnc = zeros(V_SimpleEnc.dim);
if(V_SimpleEnc.dt(1)<16)
    V_SimpleEnc.dt(1)=16; %better encoding with higher quality is needed than eg dt(1)=2. dt(1)=16 is just for safety.
end
CompEncodeOutPaths{1} = [OutputDir,filesep,'SimpleEncode.nii'];
V_SimpleEnc.fname     = CompEncodeOutPaths{1};

%Base2Encode
V_Base2Enc = V_SimpleEnc;
Y_Base2Enc = zeros(V_Base2Enc.dim);
if(V_Base2Enc.dt(1)<64)
    V_Base2Enc.dt(1)=64; %better encoding with higher quality is needed than eg dt(1)=2. dt(1)=64 is just for safety.
end
CompEncodeOutPaths{2} = [OutputDir,filesep,'Base2Encode.nii'];
V_Base2Enc.fname      = CompEncodeOutPaths{2};

%Log2Base2Encode
V_Log2Base2Encode = V_Base2Enc;
Y_Log2Base2Encode = zeros(V_Log2Base2Encode.dim);
if(V_Log2Base2Encode.dt(1)<64)
    V_Log2Base2Encode.dt(1)=64; %better encoding with higher quality
end
CompEncodeOutPaths{3}   = [OutputDir,filesep,'Log2Base2Encode.nii'];
V_Log2Base2Encode.fname = CompEncodeOutPaths{3};

%SumMap
V_SumMap = V_SimpleEnc;
Y_SumMap = zeros(V_SumMap.dim);
if(V_SumMap.dt(1)<16)
    V_SumMap.dt(1)=16; %better encoding with higher quality is needed than eg dt(1)=2. dt(1)=16 is just for safety.
end
CompEncodeOutPaths{4} = [OutputDir,filesep,'SumMap.nii'];
V_SumMap.fname        = CompEncodeOutPaths{4};

%WinnerTakeAllEncode
V_WinnerTakeAllEncode = V_SimpleEnc;
Y_WinnerTakeAllEncode = zeros(V_WinnerTakeAllEncode.dim);
if(V_WinnerTakeAllEncode.dt(1)<16)
    V_WinnerTakeAllEncode.dt(1)=16; %better encoding with higher quality is needed than eg dt(1)=2. dt(1)=16 is just for safety.
end
CompEncodeOutPaths{5}       = [OutputDir,filesep,'WinnerTakeAllEncode.nii'];
V_WinnerTakeAllEncode.fname = CompEncodeOutPaths{5};

%OverlapOnly
V_OverlapOnly = V_SimpleEnc;
Y_OverlapOnly = zeros(V_OverlapOnly.dim);
CompEncodeOutPaths{6} = [OutputDir,filesep,'OverlapOnly.nii'];
V_OverlapOnly.fname   = CompEncodeOutPaths{6};

%UniqueEncode
V_UniqueOnly = V_SimpleEnc;
Y_UniqueOnly = zeros(V_UniqueOnly.dim);
CompEncodeOutPaths{7} = [OutputDir,filesep,'UniqueOnly.nii'];
V_UniqueOnly.fname   = CompEncodeOutPaths{7};

%% processing
%Base2Encode & thresholds
Thresholds = zeros(length(ICnums),2); %pos & neg that are actually in the maps
Data_WinnerSelection = zeros(length(Y_WinnerTakeAllEncode(:)),length(ICnums)); %store all stats vals in here
for IndIC = 1:length(ICnums)
    V_tmp = spm_vol(InputFiles{IndIC});
    Y_tmp = V_tmp.private.dat(:,:,:,V_tmp.n(1));
    %get thresholds from data before applying the threshold for selection of voxels
    if(~isempty(min(Y_tmp(Y_tmp(:)>0)))) %pos threshold
        Thresholds(IndIC,1) = min(Y_tmp(Y_tmp(:)>0));
    else
        Thresholds(IndIC,1) = Inf('double');
    end
    if(~isempty(max(Y_tmp(Y_tmp(:)<0)))) %neg threshold
        Thresholds(IndIC,2) = max(Y_tmp(Y_tmp(:)<0));
    else
        Thresholds(IndIC,2) = -Inf('double');
    end
    %% threshold data
    switch(ThreshType)
        case 'P'
            Data_WinnerSelection(:,IndIC) = Y_tmp(:).*double((double(Y_tmp(:)>(max(Y_tmp(:)).*(1-ThreshP)))+double(Y_tmp(:)<(min(Y_tmp(:)).*(1-ThreshP))))>0);
        case 'Stats'
            Data_WinnerSelection(:,IndIC) = Y_tmp(:).*double((double(Y_tmp(:)>ThreshStats(1))+double(Y_tmp(:)<ThreshStats(2)))>0);
        case 'HalfMax'
            TmpY                          = Y_tmp(:).*double((double(Y_tmp(:)>ThreshStats(1))+double(Y_tmp(:)<ThreshStats(2)))>0);
            Data_WinnerSelection(:,IndIC) =  TmpY(:).*double(TmpY(:)>(max(TmpY(:))/2));
        otherwise
            error(['Unknown ThreshType: "',ThreshType,'"?!']);
    end
    
    Y_SumMap        = Y_SumMap(:)+(Data_WinnerSelection(:,IndIC)~=0); %add all significant ones
    
    Y_tmp(Data_WinnerSelection(:,IndIC)~=0) = 2^(IndIC-1); %raise nonzero parts (transformed to logicals) to current power of 2 --> base 2 encode (with next step); NB: need to start at zero.
    Y_Base2Enc      = Y_Base2Enc+Y_tmp; %adding performs the base 2 encode.
end
if((size(Y_SumMap,1)~=V_SumMap.dim(1))||(size(Y_SumMap,2)~=V_SumMap.dim(2))||(size(Y_SumMap,3)~=V_SumMap.dim(3)))
    if(length(Y_SumMap(:))~=prod(V_SumMap.dim(:)))
        error('Y_SumMap has changed size during calculations --> this is an error! Check your code & assumptions.');
    else
        Y_SumMap = reshape(Y_SumMap,V_SumMap.dim);
    end
end

%Log2Base2Encode
Y_Log2Base2Encode = Y_Base2Enc;
Y_Log2Base2Encode(Y_Log2Base2Encode>0) = log2(Y_Log2Base2Encode(Y_Log2Base2Encode>0))+1;

%OverlapOnly
Y_OverlapOnly = Y_Base2Enc;
for IndIC = 1:length(ICnums)
    Y_OverlapOnly(Y_OverlapOnly==2.^(IndIC-1)) = 0; %remove unique ones
end
Y_OverlapOnly = Y_OverlapOnly>0; %make [0,1]-mask

%UniqueOnly
Y_UniqueOnly = Y_Base2Enc;
Y_UniqueOnly(Y_OverlapOnly==1) = 0;
Y_UniqueOnly(Y_UniqueOnly>0)   = log2(Y_UniqueOnly(Y_UniqueOnly>0))+1; %each unique voxels of each component labeled by the component number

%SimpleEncode
Y_SimpleEnc = Y_Base2Enc;
Y_SimpleEnc(Y_OverlapOnly==1) = -1; %flip to negative
Y_SimpleEnc(Y_SimpleEnc>0)    = log2(Y_SimpleEnc(Y_SimpleEnc>0))+2; %each unique voxels of each component labeled by the component number+1
Y_SimpleEnc(Y_OverlapOnly==1) = 1; %bring back

%WinnerTakeAllEncode: go over nonzero voxels from SumMap
disp('creating (simple) winner take all map...');
VoxOfInterest = find(Y_SumMap>0);
for IndVox = 1:length(VoxOfInterest)
    [maximum,Y_WinnerTakeAllEncode(VoxOfInterest(IndVox))] = max(Data_WinnerSelection(VoxOfInterest(IndVox),:)); clear maximum
end    

%% write out results
%SimpleEncode
V_SimpleEnc           = spm_write_vol(V_SimpleEnc,      Y_SimpleEnc);

%Base2Encode
V_Base2Enc            = spm_write_vol(V_Base2Enc,       Y_Base2Enc);

%Log2Base2Encode
V_Log2Base2Encode     = spm_write_vol(V_Log2Base2Encode,Y_Log2Base2Encode);

%SumMap
V_SumMap              = spm_write_vol(V_SumMap         ,Y_SumMap);

%WinnerTakeAllEncode
V_WinnerTakeAllEncode = spm_write_vol(V_WinnerTakeAllEncode,Y_WinnerTakeAllEncode);

%OverlapOnly
V_OverlapOnly         = spm_write_vol(V_OverlapOnly,    Y_OverlapOnly);

%UniqueOnly
V_UniqueOnly          = spm_write_vol(V_UniqueOnly,     Y_UniqueOnly);

%Thresholds, Thresh, InputFiles
save([OutputDir,filesep,'Inputs_Thresholds_Vols.mat'],'Thresh','Thresholds','InputFiles','V_SimpleEnc','V_Base2Enc','V_Log2Base2Encode','V_SumMap','V_WinnerTakeAllEncode','V_OverlapOnly','V_UniqueOnly');

%% Done with component encode map creation.
disp('Done with creating component encode NIFTI-files.');

%% create overlap matrix
OverlapMatStruct = CreateOverlapMat(InputFiles,ThreshType,Thresh); %ie take all, do thresholding before.
save([OutputDir,filesep,'OverlapMat.mat'],'OverlapMatStruct');

%% assign DataStruct
DataStruct.InputFiles = InputFiles;
DataStruct.ICnums     = ICnums;
DataStruct.Thresh     = Thresh;
save([OutputDir,filesep,'DataStruct.mat'],'DataStruct');

%% Done.
disp('DONE.');
disp(' ');

end