function [H,SlObj,params,SliceIndices,colors,bg] = DisplayClusters(SliceSelType,NIFTI_files,xslices,colors)
%DisplayClusters.m
%
% This Script can be used to produce Slice Overlay plots from NIFTI-files.
% THIS SCRIPT IS MAINLY USEFUL FOR PLOTTING ATLASES OR MASKS.
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
%SliceSelType can either be
%    1. 'Simple'  --> Just pick all unique slices from lowest to highest slice
%                     in all nifit files that are not empty, i.e. contain values
%                     other than zero (Background). 
%                     If more then MaxNumSlices, then thin out by keeping
%                     only those that are more than minimal resolution apart.
%    2. 'Cluster' --> ASSUME that all input images are Cluster maps or masks.
%                     Pick slices that cover each cluster in each image as 
%                     best as possible, i.e. top & bottom and up to three slices
%                     with most voxels in them (sorted descending).
%                     NB: if only one cluster then we assume that is a [0,1] mask
%                     and take N slices, 
%                     top M and just above median + median + just below median + two below median if possible
%                               just above 3rdQrt + 3rdQrt + just below 3rdQrt if possible
%                               just above 1stQrt + 1stQrt + just below 1stQrt if possible 
%                     if possible.
%                     Then thin out till MaxNumSlices is reached. 
%
%BUT IT ONLY MAKES SENSE TO USE 'Cluster' HERE.
%(See SuggestSlices.m for more help.)
%
%
%USAGE:
%       [H,SlObj,params] = DisplayClusters(SliceSelType,MaskFilesPath); 
%                                          %MaskFilesPath can be empty input or not input at all, then spm_select is used. If more then one input then all are put together in one map. 
%                                          %SliceSelType can be empty or not input then the default 'Cluster' is used.
%                                          %xslices controls the number of slices per row. (can be unset or empty then a "optimum" will be determined automatically.)
%
%
%       [H,SlObj,params] = DisplayClusters(); %select slices using default option 'Cluster' %user will be asked to select nifti containing mask or clusters indicated by integers for display.
%       [H,SlObj,params] = DisplayClusters('Cluster'); %user will be asked to select nifti containing mask or clusters indicated by integers for display.
%       [H,SlObj,params] = DisplayClusters('Cluster',NIImaskPath); %display mask nifti stored at location indicated in "NIImaskPath", either as a string or cellstring und slice selection option 'Cluster'
%       [H,SlObj,params] = DisplayClusters('Simple',NIImaskPath); %display mask nifti stored at location indicated in "NIImaskPath", either as a string or cellstring und slice selection option 'Simple'
%       [H,SlObj,params] = DisplayClusters([-42,-23,0,10,20],NIImaskPath); %display slices at -42mm -23mm 0mm 10mm and 20mm in z-direction for mask nifti stored at location indicated in "NIImaskPath", either as a string or cellstring.
%       [H,SlObj,params] = DisplayClusters('Cluster',[],[],MyColors);  %Display clusters from files that are selected manually and number of slices per row is optimized and use MyColors for display.
%       
%
%
%
%V2.21
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V2.21: 01.12.2016: color assignment can be problematic when the range is not as the colors. FROM NOW ON WE ASSUME 1:NColors!!! V2.2B: 28.10.2016: changed input mode to be more robust than varargin. V2.1B: 06.10.2016: edit for special overlap display that allows input of colors instead of creation with dustinguishable_colors. V2.1: 24.02.2015: multiple inputs become mask overlays giving each input one color. V2.0: 23.04.2015: Output object and parameters for plotting and saving of plot. V1.6: 15.03.2015: Additional support of slice selection beyond 'Simple'. V1.5: (31.01.2015): Additional outputs and inputs for more control. V1.0: (14.12.2014): initial implementation (changed version from DisplayMasksNIFTIsOnSlices.m)

%% defaults
WaitTime            = 1; %s DEFAULT time for "Color" message box staying on screen if not clicked by user
Orientation         = 'axial'; %orientation of slices
TransparencyOverlay = 1;  %1==   NO   transparency NB:maybe this should be adjusted for special cases when we might have multiple overlappig maps???
structural_img      = [fileparts(which('spm.m')),filesep,'canonical',filesep,'single_subj_T1.nii'];  % standard SPM Structural or select image?    

%% Check inputs
%% select image or images for display
if(~exist('NIFTI_files','var'))
    NIFTI_files = cellstr(spm_select([1 Inf],'image','Select NIFTI-file(s) for creation of Overlay ...'));
else
    if(isempty(NIFTI_files))
        NIFTI_files = cellstr(spm_select([1 Inf],'image','Select NIFTI-file(s) for creation of Overlay ...'));
    else
        if(~iscellstr(NIFTI_files))
            if(ischar(NIFTI_files))
                NIFTI_files = cellstr(NIFTI_files);
            else
                error('NIFTI_files must be a cellstring!');
            end
        end
    end
end

%% xslices %set number of slices in a row
if(~exist('xslices','var'))
    xslices = []; %set number of slices in a row to "optimum"
else
    if(~isnumeric(xslices))
        error('xslices must be a numeric type or empty!');
    end
end

%% colors
if(~exist('colors','var'))
    colors = []; %set colors later.
else
    if(size(colors,2)~=3)
        error('colors must be a matrix indicating the colors for cluster plotting as NColors-x-3 RGB values.');
    end
end

%% check basic inputs --> slice selection type
try
    if(~ischar(SliceSelType)&&isnumeric(SliceSelType))
        SliceIndices = SliceSelType;
        disp('Slices are defined by user input.');
    else
        disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
        SliceIndices = SuggestSlices(NIFTI_files,SliceSelType); %'Cluster'); %'Simple'); %
    end
catch CATCH_SliceSelType
    disp_catch(CATCH_SliceSelType,mfilename,'CATCH_SliceSelType');
    SliceSelType = 'Cluster'; %'Simple'
    disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
    SliceIndices = SuggestSlices(NIFTI_files,SliceSelType); %'Cluster'); %'Simple'); %
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
NInputs   = length(NIFTI_files);
if(NInputs>1) %--> one color per input
    NColors = NInputs;
else
    NIItmp = nifti(NIFTI_files{1});
    if(length(NIItmp.dat.dim)==4) %4D
        [BaseTmp,FNameTmp,ExtTmp] = fileparts(NIFTI_files{1}); clear BaseTmp FNameTmp
        StartInd = regexp(ExtTmp,'\d');
        if(isempty(StartInd))
            DataExtNum = 1;
        else
            DataExtNum = str2num(ExtTmp(StartInd:end));
        end
        UniqueNonZeroVals = unique(NIItmp.dat(:,:,:,DataExtNum));
        UniqueNonZeroVals(UniqueNonZeroVals==0) = []; %remove zero (of course)
        NColors = length(UniqueNonZeroVals);
    else
        UniqueNonZeroVals = unique(NIItmp.dat(:));
        UniqueNonZeroVals(UniqueNonZeroVals==0) = []; %remove zero (of course)
        NColors = length(UniqueNonZeroVals);
    end
end

%% beautiful colors needed#111111!!!
% colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0; .5 .5 .5]; %white black gray
%bg = [1 1 1; 0 0 0; .1 .1 .1; .25 .25 .25; .5 .5 .5; .75 .75 .75; .9 .9 .9]; %white black gray(s)
%just making sure that the gray and black tones are ignored by distinguishable_colors! ;)
bg = [1 1 1; 0 0 0; .05 .05 .05; .1 .1 .1; .15 .15 .15; .25 .25 .25; .35 .35 .35; .45 .45 .45; .5 .5 .5; .65 .65 .65; .75 .75 .75; .85 .85 .85; .9 .9 .9]; %white black gray(s)

if(isempty(colors))
    %make colors
    try
        colors = distinguishable_colors(NColors,bg);
    catch CATCH_distinguishable_colors
        disp_catch(CATCH_distinguishable_colors,mfilename,'CATCH_distinguishable_colors');
        disp(['NColors= ',num2str(NColors),'.']);
        return;
    end
else
    if(size(colors,1)~=NColors)
        disp(['WARNING: Number of input colors (',num2str(size(colors,1)),') does not match the values in the input data(',num2str(NColors),').']);
    end
end
Example= repmat([size(colors,1):-1:1]',1,10); %example for plotting
figure(81); imagesc(Example); colormap(colors); title(['All ',num2str(NColors),' possible colors for plotting']); colorbar; axis('off');
TxtSize = 12; %size of the text
Order = max([1 ceil(log10(size(colors,1)))]); %number of zeros to add.
for Ind=1:size(colors,1) %NB color have to be assigned in inverse order because Example is changed.
    TxtColor = [0 0 0];
    while(all(TxtColor==colors(Ind,:))||(sqrt(sum((TxtColor-colors(Ind,:)).^2))<.15)||(sqrt(sum((TxtColor-colors(Ind,:)).^2))<.25)) %if not different enough
        TxtColor = [1 1 1]-colors(Ind,:);
        if(sqrt(sum((TxtColor-colors(Ind,:)).^2))<.15)
            TxtColor = colors(Ind,randperm(size(colors,2))); %random color change
        else
            if(sqrt(sum((TxtColor-colors(Ind,:)).^2))<.25)
                TxtColor = rand(1,size(colors,2)); %fully random colors
            end
        end
    end
    text(2,Ind,[num2str(Example(Ind,1),['%0',num2str(Order),'g']),'. rgb= [',num2str(colors(Example(Ind,1),1),'%#1.2g')],'FontSize',TxtSize,'Color',TxtColor) %NB color have to be assigned in inverse order because Example is changed.
    text(5,Ind,num2str(colors(Example(Ind,1),2),'%#1.2g'),'FontSize',TxtSize,'Color',TxtColor) %NB color have to be assigned in inverse order because Example is changed.
    text(7,Ind,[num2str(colors(Example(Ind,1),3),'%#1.2g'),']'],'FontSize',TxtSize,'Color',TxtColor) %NB color have to be assigned in inverse order because Example is changed.
end
H_info = helpdlg('These are the colors that have been generated to indicate the clusters.','Color-Generation Results');
if(isinf(WaitTime))
    uiwait(H_info);
else
    try
        uiwait(H_info,WaitTime);
        close(H_info);
    end
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
for IndInput = 1:NInputs
    if(NInputs>1)
        params{1}.img(1+IndInput).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "atlas"
        params{1}.img(1+IndInput).cmap = colors(IndInput,:);  %ColorMap of Overlay (here: "activation")
        params{1}.img(1+IndInput).range= [0.1,1.1]; %round([min(UniqueNonZeroVals),max(UniqueNonZeroVals)]); %range with probably all the right values
        params{1}.img(1+IndInput).func = 'i1(i1<=.1)=NaN;'; %remove zeros 
    else
        params{1}.img(1+IndInput).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "atlas"
        params{1}.img(1+IndInput).cmap = colors;  %ColorMap of Overlay (here: "activation")
        if(NColors==1)
            params{1}.img(1+IndInput).range= [0.1,1.1]; %bug fix for one value mask
        else
            params{1}.img(1+IndInput).range= [1, size(colors,1)]; %round([min(UniqueNonZeroVals),max(UniqueNonZeroVals)]); %first to last value
            params{1}.img(1+IndInput).func = 'i1(i1<1)=NaN;'; %remove zeros
        end
    end
    params{1}.img(1+IndInput).hold = 0; %nearest neighbor interpolation
    params{1}.img(1+IndInput).prop = TransparencyOverlay; %Transparency setting for Overlay Image
    if(TransparencyOverlay<1) %just in case we encounter this due to tweek... 
        params{1}.img(1+IndInput).type = 'truecolour'; %ie change colors to show overlap
    else
        params{1}.img(1+IndInput).type = 'split';      %ie replace Structural below with its Value/Color
    end
end
if(NInputs>1)
    if(NInputs>=3) %avoid showing too many
        params{1}.cbar = []; %don't show colorbars (too many!)
    else
        params{1}.cbar = 1+[1:NInputs];    %Only display Colorbar for Overlay
    end
else
    if(length(UniqueNonZeroVals)<=250) %6 different colors is still alright to use colorbar
        params{1}.cbar = 2;
    else
        params{1}.cbar = []; %don't show colorbars (too many!)
    end
end

%% make text
[BaseDir, fName, ext] = fileparts(NIFTI_files{1});
[BaseDir, TopDir1] = fileparts(BaseDir);
[BaseDir, TopDir2] = fileparts(BaseDir);

ResultsLast2Dirs = ['..',filesep,TopDir2,filesep,TopDir1];

params{1}.printfile    = strrep([fName,'Overlay'],' ','_');

text_annot   = cell(NInputs+1,1);
%text_annot{1}= ['Structural-Image: ',structural_img];
for IndInput = 1:NInputs
    [BaseDir, fName, ext] = fileparts(NIFTI_files{IndInput});
    if(length(fName)<45)
        text_annot{IndInput} = ['Overlay-Image ',num2str(IndInput),':   ..',ResultsLast2Dirs,filesep,fName,ext];
    else
        text_annot{IndInput} = ['Overlay-Image ',num2str(IndInput),':   ..',filesep,fName,ext];
    end
end
text_annot{1+NInputs} = ['Areas are assigned by ',num2str(NColors),' distinguishable colors.'];

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
    pause(1); %there is a strange bug that I can not nail down that disappears when waiting for a short time at this point... ??? Why??? But it works now.
    
    %% Call slover to construct SlObj for paint
    SlObj{1} = slover(params{1});
    SlObj{1}.printstr = [SlObj{1}.printstr,' -append'];

    
    %% paint slices
    paint(SlObj{1});
    drawnow;
    
    %% write annotations from SPM.mat & xSPM-Struct
    Position= [0.005,0.96,0.99,0.05]; %[0.005,0.95,0.99,0.05]
    % create annotations
    axes('Position',Position,'Visible','off');
    text(0,0.0,text_annot,'FontSize',10);
    
    clear Position
    pause(0.5);
    drawnow;
end

%% output figure handle 
H = params{1}.figure; %assign handle

%% done
disp(' ');
disp('Done.');
disp(' ');
end


