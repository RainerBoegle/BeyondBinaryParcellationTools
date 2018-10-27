function [CollectionThreshStruct] = CollectThreshSaveAs4D(ThreshMapPaths,Thresh,OutPath4DNIFTI)
% This function loads IC maps, thresholds them and writes the results to a 4D-NIFTI file.
% The settings are saved in a struct called InfoStruct.
%
% If no output path for the 4D-NIFTI file is specified, a path will be generated
% by using the directory from which all input files come from and the name of
% the output file will be "ThreshICs.nii" in that directory.
% NB: THIS WILL ONLY WORK IF THERE IS ONLY ONE DIRECTORY THAT ALL INPUTS COME FROM!
%     If this is not the case, please input also an output path.
%
%Usage:
%       CollectionThreshStruct = CollectThreshSaveAs4D(ThreshMapPaths,Thresh,OutPath4DNIFTI);
%       CollectionThreshStruct = CollectThreshSaveAs4D(ThreshMapPaths); %get everything non-zero AND generate output path for 4D-NIFTI.
%       CollectionThreshStruct = CollectThreshSaveAs4D(); %As above, but select files manually using spm_select.
%       CollectionThreshStruct = CollectThreshSaveAs4D([],[],OutPath4DNIFTI); %select files manually using spm_select AND get everything non-zero, then put it at the path specified in OutPath4DNIFTI.
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (01.11.2016): initial implementation based on previous tests

%% check inputs
%ThreshMapPaths
if(~exist('ThreshMapPaths','var'))
    ThreshMapPaths = cellstr(spm_select(inf,'image','Select IC maps NIFTI-files...'));
else
    if(isempty(ThreshMapPaths))
        ThreshMapPaths = cellstr(spm_select(inf,'image','Select IC maps NIFTI-files...'));
    else
        if(~iscellstr(ThreshMapPaths))
            if(ischar(ThreshMapPaths))
                ThreshMapPaths = cellstr(ThreshMapPaths);
            else
                error('ThreshMapPaths has to be a cellstring or convertible char-array, indicating the paths of the NIFTI files!');
            end
        end
    end
end

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

%OutPath4DNIFTI
if(~exist('OutPath4DNIFTI','var')) %need to produce it via other paths
    OutPath4DNIFTI = unique(cellfun(@fileparts,OutPath4DNIFTI,'UniformOutput',false));
    if(length(OutPath4DNIFTI)~=1)
        error('OutPath4DNIFTI automatic generation does not work if multiple paths are possible (as extracted from the input data)! Please specify the output path in that case.');
    else
        OutPath4DNIFTI = [OutPath4DNIFTI{1},filesep,'ThreshICs.nii'];
    end
else
    if(isempty(OutPath4DNIFTI))
        OutPath4DNIFTI = unique(cellfun(@fileparts,OutPath4DNIFTI,'UniformOutput',false));
        if(length(OutPath4DNIFTI)~=1)
            error('OutPath4DNIFTI automatic generation does not work if multiple paths are possible (as extracted from the input data)! Please specify the output path in that case.');
        else
            OutPath4DNIFTI = [OutPath4DNIFTI{1},filesep,'ThreshICs.nii'];
        end
    else
        if(~ischar(OutPath4DNIFTI))
            if(~iscellstr(OutPath4DNIFTI))
                error('OutPath4DNIFTI must be a string indicating the output path for the 4D-NIFTI file!');
            else
                OutPath4DNIFTI = OutPath4DNIFTI{1};
            end
        end
    end
end

%% get stats vals and voxel coordinates
Vols = spm_vol(ThreshMapPaths);
[StatsVals,XYZvox] = GetDataFromNII(Vols,Thresh);

%% create output volume 
Vout = rmfield(Vols{1},'private'); %prepare by removing any NIFTI attributes that could cause trouble.
if(Vout.dt(1)<16)
    Vout.dt(1) = 16; %make sure this is encoded at least as a 32bit float.
end
Vout.fname = OutPath4DNIFTI;

%% check if output dir exists otherwise create it
[OutDir,OutFName,OutExt] = fileparts(OutPath4DNIFTI);
if(~exist(OutDir,'dir'))
    disp(['Directory "',OutDir,'" does not exist, will create it...']);
    mkdir(OutDir);
end

%% write out data
disp(['Writing out "',OutFName,OutExt,'" to directory "',OutDir,'"...']);
reverseStr = ''; %init
for IndVol = 1:length(StatsVals)
    Y = zeros(Vout.dim(1),Vout.dim(2),Vout.dim(3));
    CurrStats = StatsVals{IndVol};
    CurrXYZvox= XYZvox{IndVol};
    for IndVox = 1:length(CurrStats)
        Y(CurrXYZvox(IndVox,1),CurrXYZvox(IndVox,2),CurrXYZvox(IndVox,3)) = CurrStats(IndVox);
    end
    Vout.n(1) = IndVol;
    spm_write_vol(Vout,Y);
    
    msg = sprintf(['Writing out Volume %0',num2str(ceil(log10(length(StatsVals)))),'d of %0',num2str(ceil(log10(length(StatsVals)))),'d...'], IndVol, length(StatsVals));
    fprintf([reverseStr, msg]); %delete the last message and print current message.
    reverseStr = repmat(sprintf('\b'), 1, length(msg)); %setup the next reverse string
end
msg = sprintf(['Writing out Volume %0',num2str(ceil(log10(length(StatsVals)))),'d of %0',num2str(ceil(log10(length(StatsVals)))),'d done.'], IndVol, length(StatsVals));
fprintf([reverseStr, msg]); %delete the last message and print current message.
disp(' ');

%% write out InfoStruct
disp(['Saving CollectionThreshStruct.mat to directory "',OutDir,'"...']);
CollectionThreshStruct = struct('ThreshMapPaths',{ThreshMapPaths},'Thresh',Thresh,'OutPath4DNIFTI',OutPath4DNIFTI);
save([OutDir,filesep,'CollectionThreshStruct.mat'],'CollectionThreshStruct');

%% Done.
disp('Done.')
disp(' ');

end