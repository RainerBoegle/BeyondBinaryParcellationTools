function [H,SlObj,params,SliceIndices] = DisplayStats(SliceSelType,NIFTI_files,Thresholds,xslices,AddText,UseTransparency,UseColorBars)
%DisplayStats.m
%
% This Script can be used to produce Slice Overlay plots from NIFTI-files.
% THIS SCRIPT IS MAINLY USEFUL FOR PLOTTING STATISTICS IMAGES.
%
% The SPM toolbox "slover" originally created by Matthew Brett (http://imaging.mrc-cbu.cam.ac.uk/imaging/MatthewBrett)
% is used to produce the Overlays.
%
% The main purpose of this Script is to make the use of "slover" easier.
%
% It is expected that the user inputs one or more masks, i.e. volumes with
% integer valued only, that are then displayed in as many colors as there
% are inputs or unique values of the input, in case there is only one input.
%
% All parameters are fixed (more or less).
%
% The mode of suggesting slices can be set in the following ways:
%SliceSelType Mode:
%    1. 'Simple'  --> Just pick all unique slices from lowest to highest slice
%                     in all nifit files that are not empty, i.e. contain values
%                     other than zero (Background). 
%                     If more then MaxNumSlices, then thin out by keeping
%                     only those that are more than minimal resolution apart.
%    2. 'Stats'   --> ASSUME that all input images are Statistics maps.
%                     For each image find the maximum statistics value and remove
%                     the two neighboring slices and repeat with the remaining slices.
%                     Do the same starting with the minimum and ascending.
%                     Thin out each list for all images till MaxNumSlices is reached.
%
% THE BEST CHOICE HERE IS USUALLY 'Stats'.
%
% USAGE:
%       [H,SlObj,params] = DisplayStats(SliceSelType,NIFTI_files,Thresholds,xslices,AddText,UseTransparency); 
%                                       %SliceSelType can be empty or not input then the default 'Cluster' is used.
%                                       %StatsFilesPaths can be empty or not input then spm_select is used. 
%                                       %Thresholds must input image [Negative, Positive; Negative, Positive; ...; Negative, Positive]; ... 
%                                       %xslices controls the number of slices per row. (can be unset or empty then a "optimum" will be determined automatically.)
%                                       %AddText controls if info text indicating paths of NIFTI files is printed or not. If == 1 then text is printed.
%                                       %UseTransparency controls if non-significant parts of the plot are shown in semitransparent way or if they are removed completely like zero-values.
%
%       [H,SlObj,params] = DisplayStats(); %select images & default slice selection type to 'Stats'.
%       [H,SlObj,params] = DisplayStats('Simple'); %select images & default slice selection type to 'Simple'.
%       [H,SlObj,params] = DisplayStats('Stats',StatsFilesPaths); %load images from StatsFilesPaths (either string or cellstring) & default slice selection type to 'Stats'.
%       [H,SlObj,params] = DisplayStats([-42,-23,0,10,20],StatsFilesPaths); %display slices at -42mm -23mm 0mm 10mm and 20mm in z-direction %load images from StatsFilesPaths (either string or cellstring).
%       [H,SlObj,params] = DisplayStats([-42,-23,0,10,20],StatsFilesPaths,[-1.75,5]); %as above but using thresholds for display such that parts below threshold are transparent. Display slices at -42mm -23mm 0mm 10mm and 20mm in z-direction %load images from StatsFilesPaths (either string or cellstring).
%       [H,SlObj,params] = DisplayStats([-42,-23,0,10,20],StatsFilesPaths,[-1.75,5],6); %as above using 6 slices per row in display and using thresholds for display such that parts below threshold are transparent. Display slices at -42mm -23mm 0mm 10mm and 20mm in z-direction %load images from StatsFilesPaths (either string or cellstring).
%
%
%
%
%V2.6
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V2.6 (06.12.2016): allow suppression of colorbars. V2.5 (21.11.2016): New way of controlling the inputs and also set default to NON-transparent display excluding non-significant parts, but transparent plot is possible. V2.2 (16.11.2016): allow text to be printed OR not. V2.1 (24.04.2015): interpolate the transparency until threshold in 10-steps. V2.0: (23.04.2015) Output object and parameters for plotting and saving of plot. V1.5: (16.02.2015): Additional outputs and inputs for more control. V1.0: (14.12.2014): initial implementation (changed version from DisplayMasksNIFTIsOnSlices.m)

%% defaults
Orientation         = 'axial'; %orientation of slices
%TransparencyOverlay = 0.9; %0.5;  %1==   NO   transparency 
InterpolationSetting = 0; %nearest neighbor; %1; %trilinear %2->127     Higher order Lagrange (polynomial) interpolation using different holds (second-order upwards). %-127 - -1   Different orders of sinc interpolation. 
structural_img      = [fileparts(which('spm.m')),filesep,'canonical',filesep,'single_subj_T1.nii'];  % standard SPM Structural or select image?    
PosColor            = hot(128);
PosColor            = [bsxfun(@times,repmat(PosColor(20,:),19,1),[1,3,3.5,4,(5:19)]'./20);PosColor(20:end,:)];
NegColor            = winter(128);
NegColor            = NegColor(end:-1:1,:);

%% Check inputs
if(~exist('NIFTI_files','var'))
    NIFTI_files = cellstr(spm_select([1 2],'image','Select NIFTI-file(s) for creation of Overlay ...'));
elseif(isempty(NIFTI_files))
    NIFTI_files = cellstr(spm_select([1 2],'image','Select NIFTI-file(s) for creation of Overlay ...'));
else
    if(~iscellstr(NIFTI_files))
        if(ischar(NIFTI_files))
            NIFTI_files = cellstr(NIFTI_files);
        else
            error('NIFTI_files has to be a cellstring or convertible char-array.');
        end
    end
end
if(length(NIFTI_files)>2)
    error('Only two inputs allowed!');
end

if(~exist('Thresholds','var'))
    Thresholds  = nan(length(NIFTI_files),2); %no thresholds!
elseif(isempty(Thresholds))
    Thresholds  = nan(length(NIFTI_files),2); %no thresholds!
else
    if(size(Thresholds,1)~=length(NIFTI_files))
        if(size(Thresholds,2)==1&&size(Thresholds,1)==2)
            Thresholds = Thresholds';
        else
            error(['Number of input thresholds (',num2str(size(Thresholds,1)),') do not correspond to number of input files (',num2str(length(NIFTI_files)),')!']);
        end
    end
    if(size(Thresholds,2)~=2)
        error('Only two thresholds values [negative positive] possilbe!');
    end
end

if(~exist('xslices','var'))
    xslices     = []; %set number of slices in a row to "optimum"
elseif(~isempty(xslices))
    if(mod(length(xslices),1)~=0)
        error('xslices settings must be an integer!');
    end
end

if(~exist('AddText','var'))
    AddText     = 1; %add text to figure
elseif(isempty(AddText))
    AddText     = 1; %add text to figure
else
    if(AddText~=0)
        AddText     = 1; %add text to figure
    end
end

if(~exist('UseTransparency','var'))
    UseTransparency = 0; %if == 1 then display non-significant parts transparent.
elseif(isempty(UseTransparency))
    UseTransparency = 0; %if == 1 then display non-significant parts transparent.
else
    if(UseTransparency~=0)
        UseTransparency = 1; %if == 1 then display non-significant parts transparent.
    end
end

if(~exist('UseColorBars','var'))
    UseColorBars = 1;
elseif(isempty(UseColorBars))
    UseColorBars = 1;
else
    if(UseColorBars~=0)
        UseColorBars = 1;
    end
end

%% check basic inputs --> slice selection type
try
    if(~ischar(SliceSelType)&&isnumeric(SliceSelType))
        SliceIndices = SliceSelType;
        disp('Slices are defined by user input.');
    else
        disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
        SliceIndices = SuggestSlices(NIFTI_files,SliceSelType,'Thresh',Thresholds); %'Cluster'); %'Simple'); %
    end
catch CATCH_SliceSelType
    disp_catch(CATCH_SliceSelType,mfilename,'CATCH_SliceSelType');
    SliceSelType = 'Stats'; %'Simple'
    disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
    SliceIndices = SuggestSlices(NIFTI_files,SliceSelType,'Thresh',Thresholds); %'Cluster'); %'Simple'); %
end

%% Slices selection & orientation 
SliceIndicesStr = '[';
for IndSlice = 1:length(SliceIndices)
    if(IndSlice==1)
        SliceIndicesStr = [SliceIndicesStr,num2str(SliceIndices(IndSlice))];
    else
        SliceIndicesStr = [SliceIndicesStr,',',num2str(SliceIndices(IndSlice))];
    end
end
SliceIndicesStr = [SliceIndicesStr,']'];
disp(['Will produce overlay using ',Orientation,' slices ',SliceIndicesStr,'.']);

%% check if Structural image can be reached
if(~exist(structural_img))
    disp(['Structural image (',structural_img,') not found! Check paths.']);
    pause(2);
    if((evalin('base','exist(''prevsect'')')==1))
        structural_img = evalin('base','prevsect');
    else
        structural_img = spm_select(1,'image','Select structural image for background...');
    end
    if(strcmp(structural_img((length(structural_img)-1):end),',1'))
        structural_img = structural_img(1:(length(structural_img)-2));
    end
end
disp(['Using "',structural_img,'" as structural background img']);

%% check inputs and select colors
NInputs= length(NIFTI_files);
Ranges = cell(length(NIFTI_files),2);
for IndInput = 1:NInputs
    vol = spm_vol(NIFTI_files{IndInput});
    if(length(vol)>1) %4D file but shouldn't be!
        error(['Input ',num2str(IndInput),': "',NIFTI_files{IndInput},'" is 4D, but no specific image was indicated. Use spm_select to get the image from the 4D-file that you want.']);
    end
    if(length(vol.private.dat.dim)==4)
        DataNeg = vol.private.dat(:,:,:,vol.n(1)); %this is actually all data, will separate next
        DataPos = DataNeg(find(DataNeg(:)>0))';
        DataNeg = DataNeg(find(DataNeg(:)<0))';
    else
        DataNeg = vol.private.dat(find(vol.private.dat(:)<0))';
        DataPos = vol.private.dat(find(vol.private.dat(:)>0))';
    end
    
    if(~isempty(DataNeg))
        RangeNeg = minmax(DataNeg);
    else
        RangeNeg = [];
    end
    if(~isempty(DataPos))
        RangePos = minmax(DataPos);
    else
        RangePos = [];
    end
    Ranges{IndInput,1} = RangeNeg;
    Ranges{IndInput,2} = RangePos;
end

%% Fill params for slover call: "ACTIVATION OVERLAY FROM NIFTI"
params = cell(1,1);
params{1}.slices    = SliceIndices; %Slice Indices of Slices to display
params{1}.transform = Orientation;  %Slice Orientation
if(~isempty(xslices))
    params{1}.xslices = xslices; %number of slices per row.
end

%structural image as background
params{1}.img(1).vol   = spm_vol(structural_img);  %get Structural Image
params{1}.img(1).cmap  = gray(256); %ColorMap of Structural Image
params{1}.img(1).range = minmax(params{1}.img(1).vol.private.dat(:)'); %Displayed Range of Structural Image i.e. all values
params{1}.img(1).prop  = 1; %Transparency setting for Structural Image
params{1}.img(1).type  = 'truecolour'; %Image which can be overlayed by other image.

%overlays
N = 1; %init
for IndInput = 1:NInputs
    if(~isempty(Ranges{IndInput,1})&&~isempty(Ranges{IndInput,2})) %negative & positive  
        %negative
        params{1}.img(N+1).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+1).cmap = NegColor;  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+1).range= [min(Ranges{IndInput,1}(:)),0]; %,max(Ranges{IndInput}(:))]; %range for colors
        params{1}.img(N+1).hold = InterpolationSetting; %.hold = 0; %nearest neighbor interpolation
        if(~isnan(Thresholds(IndInput,1))) %negative
            if(UseTransparency)
                for IndTransp = 1:9
                    params{1}.img(N+1).func = ['i1(i1<',num2str(Thresholds(IndInput,1).*(IndTransp/10)),'|i1>=0)=NaN;']; %remove range below negative threshold and above zeros
                    params{1}.img(N+1).prop = IndTransp/10; %Transparency setting for Overlay Image
                    params{1}.img(N+1).type = 'truecolour'; %ie change colors to show overlap
                    N = length(params{1}.img);
                    if(IndTransp<9)
                        params{1}.img(N+1)=params{1}.img(N);
                    end
                end
                %non-transparent part
                params{1}.img(N+1)=params{1}.img(N);
                params{1}.img(N+1).func = ['i1(i1>',num2str(Thresholds(IndInput,1)),')=NaN;']; %remove range above threshold (negative)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            else
                params{1}.img(N+1).func = ['i1(i1>',num2str(Thresholds(IndInput,1)),')=NaN;']; %remove range above threshold (negative)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            end
        else
            params{1}.img(N+1).func = ['i1(i1>=0)=NaN;']; %remove range above zeros
            params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
            params{1}.img(N+1).type = 'split';      %ie replace Structural below with its Value/Color
            N = length(params{1}.img);
        end
        
        %positive
        params{1}.img(N+1).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+1).cmap = PosColor;  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+1).range= [0,max(Ranges{IndInput,2}(:))]; %range for colors
        params{1}.img(N+1).hold = InterpolationSetting; %.hold = 0; %nearest neighbor interpolation
        if(~isnan(Thresholds(IndInput,2))) %positive
            if(UseTransparency)
                for IndTransp = 1:9
                    params{1}.img(N+1).func = ['i1(i1>',num2str(Thresholds(IndInput,2).*(IndTransp/10)),'|i1<=0)=NaN;']; %remove range below zeros and above threshold (positive)
                    params{1}.img(N+1).prop = IndTransp/10; %Transparency setting for Overlay Image
                    params{1}.img(N+1).type = 'truecolour'; %ie change colors to show overlap
                    N = length(params{1}.img);
                    if(IndTransp<9)
                        params{1}.img(N+1)=params{1}.img(N);
                    end
                end
                %non-transparent part
                params{1}.img(N+1)=params{1}.img(N);
                params{1}.img(N+1).func = ['i1(i1<',num2str(Thresholds(IndInput,2)),')=NaN;']; %remove range below threshold (positive)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            else
                params{1}.img(N+1).func = ['i1(i1<',num2str(Thresholds(IndInput,2)),')=NaN;']; %remove range below threshold (positive)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            end
        else
            params{1}.img(N+1).func = ['i1(i1<=0)=NaN;']; %remove range above zeros
            params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
            params{1}.img(N+1).type = 'split';      %ie replace Structural below with its Value/Color
            N = length(params{1}.img);
        end
    elseif(~isempty(Ranges{IndInput,1})) %negative
        params{1}.img(N+1).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+1).cmap = NegColor;  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+1).range= [min(Ranges{IndInput,1}(:)),0]; %,max(Ranges{IndInput}(:))]; %range for colors
        params{1}.img(N+1).hold = InterpolationSetting; %.hold = 0; %nearest neighbor interpolation
        if(~isnan(Thresholds(IndInput,1))) %negative
            if(UseTransparency)
                for IndTransp = 1:9
                    params{1}.img(N+1).func = ['i1(i1<',num2str(Thresholds(IndInput,1).*(IndTransp/10)),'|i1>=0)=NaN;']; %remove range below negative threshold and above zeros
                    params{1}.img(N+1).prop = IndTransp/10; %Transparency setting for Overlay Image
                    params{1}.img(N+1).type = 'truecolour'; %ie change colors to show overlap
                    N = length(params{1}.img);
                    if(IndTransp<9)
                        params{1}.img(N+1)=params{1}.img(N);
                    end
                end
                %non-transparent part
                params{1}.img(N+1)=params{1}.img(N);
                params{1}.img(N+1).func = ['i1(i1>',num2str(Thresholds(IndInput,1)),')=NaN;']; %remove range above threshold (negative)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            else
                params{1}.img(N+1).func = ['i1(i1>',num2str(Thresholds(IndInput,1)),')=NaN;']; %remove range above threshold (negative)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            end
        else
            params{1}.img(N+1).func = ['i1(i1>=0)=NaN;']; %remove range above zeros
            params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
            params{1}.img(N+1).type = 'split';      %ie replace Structural below with its Value/Color
            N = length(params{1}.img);
        end
    elseif(~isempty(Ranges{IndInput,2})) %positive
        params{1}.img(N+1).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+1).cmap = PosColor;  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+1).range= [0,max(Ranges{IndInput,2}(:))]; %range for colors
        params{1}.img(N+1).hold = InterpolationSetting; %.hold = 0; %nearest neighbor interpolation
        if(~isnan(Thresholds(IndInput,2))) %positive
            if(UseTransparency)
                for IndTransp = 1:9
                    params{1}.img(N+1).func = ['i1(i1>',num2str(Thresholds(IndInput,2).*(IndTransp/10)),'|i1<=0)=NaN;']; %remove range below zeros and above threshold (positive)
                    params{1}.img(N+1).prop = IndTransp/10; %Transparency setting for Overlay Image
                    params{1}.img(N+1).type = 'truecolour'; %ie change colors to show overlap
                    N = length(params{1}.img);
                    if(IndTransp<9)
                        params{1}.img(N+1)=params{1}.img(N);
                    end
                end
                %non-transparent part
                params{1}.img(N+1)=params{1}.img(N);
                params{1}.img(N+1).func = ['i1(i1<',num2str(Thresholds(IndInput,2)),')=NaN;']; %remove range below threshold (positive)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            else
                params{1}.img(N+1).func = ['i1(i1<',num2str(Thresholds(IndInput,2)),')=NaN;']; %remove range below threshold (positive)
                params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
                params{1}.img(N+1).type = 'split'; %ie change colors to show overlap
                N = length(params{1}.img);
            end
        else
            params{1}.img(N+1).func = ['i1(i1<=0)=NaN;']; %remove range above zeros
            params{1}.img(N+1).prop = 1; %Transparency setting for Overlay Image
            params{1}.img(N+1).type = 'split';      %ie replace Structural below with its Value/Color
            N = length(params{1}.img);
        end
    end
end
%colorbars?
if(UseColorBars)
    if(length(params{1}.img)==3)
        params{1}.cbar = 2:3;
    elseif(length(params{1}.img)>3)
        params{1}.cbar = [2,(length(params{1}.img)-1)];
    else
        params{1}.cbar = 2;
    end
    if(NInputs>1)
        params{1}.cbar = []; %Problem(need common colorbar!): this is a bug fix, need to think carefully about multiple inputs to make a better solution!!!
    end
else
    params{1}.cbar = []; %NO colorbar
end

%% make text
[BaseDir, fName, ext] = fileparts(NIFTI_files{1});
[BaseDir, TopDir1] = fileparts(BaseDir);
[BaseDir, TopDir2] = fileparts(BaseDir);

ResultsLast2Dirs = [filesep,TopDir2,filesep,TopDir1];
ResultsLastDir   = [filesep,TopDir1];

params{1}.printfile    = strrep([fName,'Overlay'],' ','_');

text_annot   = cell(NInputs,1);
%text_annot{1}= ['Structural-Image: ',structural_img];
for IndInput = 1:NInputs
    [BaseDir, fName, ext] = fileparts(NIFTI_files{IndInput});
    if(length(fName)>45)
        text_annot{IndInput} = ['Overlay-Image ',num2str(IndInput),':   ..',filesep,fName,ext];
    elseif(length(fName)>20)
        text_annot{IndInput} = ['Overlay-Image ',num2str(IndInput),':   ..',ResultsLastDir,filesep,fName,ext];
    else
        text_annot{IndInput} = ['Overlay-Image ',num2str(IndInput),':   ..',ResultsLast2Dirs,filesep,fName,ext];
    end
end

clear BaseDir fName TopDir1 TopDir2 ResultsLast2Dirs % Step_cmap thresInd

%clear to avoid pain on rerun. ;)
clear check Orientation structural_img 

%% make overlay & ask for saveing or not
%% create slover object & "print" to graphics window
SlObj = cell(1,1);
if(isfield(params{1},'img'))
    %% add a new window and try to leave space for the text 
    params{1}.figure        = spm_figure(); %new figure
    params{1}.area.position = [0.005,0.005,0.99,0.95];
    params{1}.area.units    = 'normalized';
    params{1}.area.halign   = 'center';
    params{1}.area.valign   = 'middle';
    pause(1);
    
    %% Call slover to construct SlObj for paint
    SlObj{1} = slover(params{1});
    if(1)%for debugging
        SlObj{1}.printstr = [SlObj{1}.printstr,' -append'];
    end
    
    %% paint slices
%     spm_figure('GetWin','Graphics');
    paint(SlObj{1});
    drawnow;
    
    %% write annotations from SPM.mat & xSPM-Struct
    if(AddText)
        Position= [0.005,0.96,0.99,0.05]; %[0.005,0.95,0.99,0.05]
        % create annotations
        axes('Position',Position,'Visible','off');
        text(0,0.0,text_annot,'FontSize',10);
        
        clear Position
        pause(0.5);
        drawnow;
    end
end

%% output figure handle 
H = params{1}.figure;


%% done
disp(' ');
disp('Done.');
disp(' ');
end



