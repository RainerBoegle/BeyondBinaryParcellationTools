function OutPath = saveThresMap(MapExtractStruct,OutPath)
% This function can be used to output a NIFTI-file from the Parameters collected 
% (from a NIFTI) using GetParamsFromMap.m OR (from a SPM.mat contrast) using GetParamsFromSPMmat.m, 
% and that were saved in the struct "MapExtractStruct" as a NIFTI file, i.e. a thresholded map.
%
%       MapExtractStruct.
%                       .Thresholds   <--   Positive & Negative threshold for Map (can be "eps" if >0 or <0 is wanted).
%
%                       .Negative.    <--   Negative threshold data (if used, otherwise empty)
%                       .Positive.    <--   Positive threshold data
%                                .Coords_mm   <<   The mm-coordinates of all significant voxels, i.e. "above" threshold.
%                                .StatsVals   <<   The statistics values of all significant voxels, i.e. "above" threshold.
%                                .Coords_vox  <<   The voxel-coordinates of all significant voxels, i.e. "above" threshold.
%
%                       .V_map.       <--   SPM-vol struct of the input map (can be used for later output of NIFTI).
%
%USAGE:
%       OutPath = saveThresMap(MapExtractStruct,OutPath);
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (27.04.2015): initial implementation

%% check inputs
try
    if(ischar(MapExtractStruct)) %filepath?
        if(exist(MapExtractStruct,'file'))
            load(MapExtractStruct,'MapExtractStruct');
        else
            error('Input has a to be a file path to the "MapExtractStruct" or the "MapExtractStruct" itself.');
        end
    else
        if(isstruct(MapExtractStruct))
            if(~isfield(MapExtractStruct,'V_map'))
                error('The input "MapExtractStruct" does not contain V_map field!');
            end
        else
            if(isempty(MapExtractStruct))
                load(spm_select(1,'mat','Select *.mat-file containing "MapExtractStruct" structure...'));
            else
                error('Input has a to be a file path to the "MapExtractStruct" or the "MapExtractStruct" itself.');
            end
        end
    end
catch CATCH_MapExtractStruct
    disp_catch(CATCH_MapExtractStruct,mfilename,'CATCH_MapExtractStruct')
    load(spm_select(1,'mat','Select *.mat-file containing "MapExtractStruct" structure...'));
end

try
    if(~isempty(fileparts(OutPath)))
        if(~exist(fileparts(OutPath),'dir'))
            error(['Directory "',fileparts(OutPath),'" DOES NOT EXIST!!!']);
        end
    else
        if(isempty(OutPath))
            [FileName,PathName,FilterIndex] = uiputfile('*.nii','Save THRESHOLDED Stats-Map?','THRES_StatsMap.nii');
            if(FilterIndex~=0)
                OutPath = [PathName,filesep,FileName];
            end
        else
            if(ischar(OutPath))
                OutPath = [pwd,filesep,OutPath];
            else
                [FileName,PathName,FilterIndex] = uiputfile('*.nii','Save THRESHOLDED Stats-Map?','THRES_StatsMap.nii');
                if(FilterIndex~=0)
                    OutPath = [PathName,filesep,FileName];
                end
            end
        end
    end
catch CATCH_OutPath
    disp_catch(CATCH_OutPath,mfilename,'CATCH_OutPath')
    [FileName,PathName,FilterIndex] = uiputfile('*.nii','Save THRESHOLDED Stats-Map?','THRES_StatsMap.nii');
    if(FilterIndex~=0)
        OutPath = [PathName,filesep,FileName];
    end
end

%% make new volume information
V_ThresMapOut = MapExtractStruct.V_map;
V_ThresMapOut.fname = OutPath; %assign output filename and path

Map_dat = zeros(MapExtractStruct.V_map.dim); %init empty ie zeros
if(~isempty(MapExtractStruct.Negative)) %fill with negative stats values if available
    for IndVox = 1:size(MapExtractStruct.Negative.Coords_vox,1)
        Map_dat(MapExtractStruct.Negative.Coords_vox(IndVox,1),MapExtractStruct.Negative.Coords_vox(IndVox,2),MapExtractStruct.Negative.Coords_vox(IndVox,3)) = MapExtractStruct.Negative.StatsVals(IndVox);
    end
end
if(~isempty(MapExtractStruct.Positive)) %fill with positive stats values if available
    for IndVox = 1:size(MapExtractStruct.Positive.Coords_vox,1)
        Map_dat(MapExtractStruct.Positive.Coords_vox(IndVox,1),MapExtractStruct.Positive.Coords_vox(IndVox,2),MapExtractStruct.Positive.Coords_vox(IndVox,3)) = MapExtractStruct.Positive.StatsVals(IndVox);
    end
end

if(V_ThresMapOut.dt(1)<16)
    V_ThresMapOut.dt(1) = 16; %not necessary but save
end
if(V_ThresMapOut.n(1)~=1) %if map is 4D
    V_ThresMapOut.n(1) = 1; %3D file has only index 1 as forth nothing else.
end
if(isfield(V_ThresMapOut,'private'))
    if(length(V_ThresMapOut.private.dat.dim)==4) %if map is 4D originally then we need to remove this information otherwise spm_write_vol will try to write this as 4D
        V_ThresMapOut = rmfield(V_ThresMapOut,'private');
    end
end

%% write NIFTI
V_ThresMapOut = spm_write_vol(V_ThresMapOut,Map_dat);

%% Done.
[OutDir,OutfName,Outext] = fileparts(V_ThresMapOut.fname);
disp(' ');
disp(['Thresholded Statistics-Map has been written to NIFTI-file "',OutfName,Outext,'".']);
disp(['In the directory "',OutDir,'".']);

end