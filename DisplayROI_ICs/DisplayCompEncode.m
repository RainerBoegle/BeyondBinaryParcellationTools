function [H,SlObj,params,SliceIndices,AllColors,AllRanges] = DisplayCompEncode(SliceSelType,CompEncodeDir,xslices,DispColors,FuncVerLog2B2E_ColGen,RenderOverlays)
%DisplayCompEncode.m
%
% This Script can be used to produce Slice Overlay plots from CompEncode NIFTI-files.
% THIS SCRIPT IS MAINLY USEFUL FOR PLOTTING ATLASES OR MASKS THAT ARE IN CompEncode.
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
%       [H,SlObj,params,SliceIndices,AllColors,AllRanges] = DisplayCompEncode(SliceSelType,CompEncodeDir,xslices); 
%                                                           %SliceSelType can be empty or not input then the default 'Cluster' is used.
%                                                           %CompEncodeDir can be empty input or not input at all, then spm_select is used. This specifies the directory with the ComponentEncode NIFTI-files. 
%                                                           %xslices controls the number of slices per row. (can be unset or empty then a "optimum" will be determined automatically.)
%
%
%       [H,SlObj,params] = DisplayCompEncode(); %select slices using default option 'Cluster' %user will be asked to select nifti containing mask or clusters indicated by integers for display.
%       [H,SlObj,params] = DisplayCompEncode('Cluster'); %user will be asked to select nifti containing mask or clusters indicated by integers for display.
%       [H,SlObj,params] = DisplayCompEncode('Cluster',CompEncodeDir); %display niftis stored in directory indicated in "CompEncodeDir", -a string; 2.slice selection option 'Cluster'
%       [H,SlObj,params] = DisplayCompEncode('Simple' ,CompEncodeDir); %display niftis stored in directory indicated in "CompEncodeDir", -a string; 2.slice selection option 'Simple'
%       [H,SlObj,params] = DisplayCompEncode([-42,-23,0,10,20],CompEncodeDir); %display slices at -42mm -23mm 0mm 10mm and 20mm in z-direction for niftis stored in directory indicated in "CompEncodeDir", -a string.
%       [H,SlObj,params] = DisplayCompEncode('Cluster',CompEncodeDir,6); %show 6 slices per row (fixed) and don't optimize display! 
%       
%
%
%
%V1.1
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.1: (26.10.2016): add option to automatically render the overlays V1.0: (25.08.2015): initial implementation (changed version from DisplayClusters.m)

%% defaults
Orientation    = 'axial'; %orientation of slices
OpacityOverlay = 1;  %1==   NO   transparency NB:maybe this should be adjusted for special cases when we might have multiple overlappig maps???
structural_img = [fileparts(which('spm.m')),filesep,'canonical',filesep,'single_subj_T1.nii'];  % standard SPM Structural or select image?    

%% check inputs
%CompEncodeDir --> get nifti-files
try
    if(~exist(CompEncodeDir))
        CompEncodeDir = spm_select(1,'dir','Select directory with component encode NIFTI-files...');
    else
        %search for files
        NIFTI_files = GetCompEncodeFiles(CompEncodeDir);
    end
catch CATCH_CompEncodeDir
    disp_catch(CATCH_CompEncodeDir,[mfilename,'>CATCH_CompEncodeDir'],'CATCH_CompEncodeDir'); %assign & display error as warning
    CompEncodeDir = spm_select(1,'dir','Select directory with component encode NIFTI-files...');
    %search for files
    NIFTI_files = GetCompEncodeFiles(CompEncodeDir);
end
NInputs = length(NIFTI_files);
%SliceSelType
try
    if(~ischar(SliceSelType)&&isnumeric(SliceSelType))
        SliceIndices = SliceSelType;
        disp('Slices are defined by user input.');
    else
        disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
        SliceIndices = SuggestSlices(NIFTI_files,SliceSelType); %'Cluster'); %'Simple'); %
    end
catch CATCH_SliceSelType
    disp_catch(CATCH_SliceSelType,[mfilename,'>CATCH_SliceSelType'],'CATCH_SliceSelType');
    SliceSelType = 'Cluster'; %'Simple'
    disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
    SliceIndices = SuggestSlices(NIFTI_files,SliceSelType); %'Cluster'); %'Simple'); %
end
%xslices   = [];
try
    if(isempty(xslices))
        disp('Will optimize display of slices to fill figure.');
    else
        if(~isvector(xslices)||~isnumeric(xslices))
            error('xslices must be a numeric vector indicating the slices in mm that should be displayed.');
        else
            disp(['Will use ',num2str(xslices),' slices per row for display.']);
        end
    end
catch CATCH_xslices
    disp_catch(CATCH_xslices,[mfilename,'>CATCH_xslices'],'CATCH_xslices'); %assign & display error as warning
    xslices   = [];
    disp('Will optimize display of slices to fill figure. (default)');
end
%DispColors = 1; %display colors in extra figure
try
    if(isempty(DispColors))
        disp('Not displaying colors separately, good luck figuring out which is which...');
        DispColors = 0;
    else
        if(DispColors==0)
            disp('Not displaying colors separately, good luck trying to see which is which...');
        end
    end
catch CATCH_DispColors
    disp_catch(CATCH_DispColors,[mfilename,'>CATCH_DispColors'],'CATCH_DispColors'); %assign & display error as warning
    DispColors = 1; %Will display colors per default.
end  
try
    FuncVerLog2B2E_ColGen;
catch CATCH_FuncVerLog2B2E_ColGen
    disp_catch(CATCH_FuncVerLog2B2E_ColGen,[mfilename,'>CATCH_FuncVerLog2B2E_ColGen'],'CATCH_FuncVerLog2B2E_ColGen'); %assign & display error as warning
    FuncVerLog2B2E_ColGen = 1;
end

%RenderOverlays
if(~exist('RenderOverlays','var'))
    RenderOverlays = nan;
end
if(~isnan(RenderOverlays))
    if(RenderOverlays)
        disp('Will automatically render all overlays...');
    end
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

%% Fill params for slover call: "OVERLAY FROM NIFTI"
%  One overlay per input NIFTI-file.
%  Therefore also need one set per input NIFTI-file.

disp('Creating colors and settings for overlays...');
H        = cell(NInputs,1);
SlObj    = cell(NInputs,1);
params   = cell(NInputs,1);
AllColors= cell(NInputs,1);
AllRanges= cell(NInputs,1);
NamesForRendering = cell(NInputs,1);
for IndInput = 1:NInputs
    [tmp,ColorGenType] = fileparts(NIFTI_files{IndInput}); clear tmp
    switch(ColorGenType)
        case 'Base2Encode'
            disp('Skipping over Base2Encode...');
            continue;
        case 'SumMap'
            disp(ColorGenType);
            %special case for SumMap --> produce histogram
            PlotSumMapHist(NIFTI_files{IndInput});
        otherwise
            disp(ColorGenType);
    end
    NamesForRendering{IndInput} = [ColorGenType,'_',datestr(now,'ddmmmyyyy_HHMM')];
    
    %slices & orientation
    params{IndInput}.slices    = SliceIndices; %Slice Indices of Slices to display
    params{IndInput}.transform = Orientation;  %Slice Orientation
    if(~isempty(xslices))
        params{IndInput}.xslices = xslices; %number of slices per row.
    end
    
    %structural image as background
    params{IndInput}.img(1).vol   = spm_vol(structural_img);  %get Structural Image
    params{IndInput}.img(1).cmap  = gray(256); %ColorMap of Structural Image
    params{IndInput}.img(1).range = minmax(params{IndInput}.img(1).vol.private.dat(:)'); %Displayed Range of Structural Image i.e. all values
    params{IndInput}.img(1).prop  = OpacityOverlay; %Transparency setting for Structural Image
    params{IndInput}.img(1).type  = 'split'; %'truecolour'; %Image which can be overlayed by other image.
    
    %create colors & range
    [Colors,ColorRange,ColorInfoStr,NComp] = CreateColors(NIFTI_files{IndInput},ColorGenType,DispColors,FuncVerLog2B2E_ColGen);
    AllColors{IndInput} = Colors;     %assign for output
    AllRanges{IndInput} = ColorRange; %assign for output
    
    %overlays
    params{IndInput}.img(2).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "atlas"
    params{IndInput}.img(2).cmap = Colors;  %ColorMap of Overlay (here: "activation")
    params{IndInput}.img(2).range= ColorRange; %round([min(UniqueNonZeroVals),max(UniqueNonZeroVals)]); %range with probably all the right values
    params{IndInput}.img(2).func = 'i1(i1<=0.1)=NaN;'; %remove zeros
    switch(ColorGenType)
        case 'Log2Base2Encode'
            if(NComp<40)
                params{IndInput}.img(2).hold  = -3; %1 is trilinear; >1 Lagrange polynomial interpolation order; -1...-127 sinc interpolation order
            else
                if(NComp<60)
                    params{IndInput}.img(2).hold  = 1; %-3; %1 is trilinear; >1 Lagrange polynomial interpolation order; -1...-127 sinc interpolation order
                else
                    params{IndInput}.img(2).hold  = 0; %-3; %1 is trilinear; >1 Lagrange polynomial interpolation order; -1...-127 sinc interpolation order
                end
            end
            params{IndInput}.img(2).type  = 'split'; %
        otherwise
            params{IndInput}.img(2).hold  = 0; %nearest neighbor interpolation
            params{IndInput}.img(2).type  = 'split'; %
    end
    
    %colorbar
    params{IndInput}.cbar = 2; %2==show colorbar of CompEncode %[]==don't show colorbars 
    
    
    %% make text
    [BaseDir, fName, ext] = fileparts(NIFTI_files{IndInput});
    [tmpDir, TopDir1] = fileparts(BaseDir);
    [tmpDir, TopDir2] = fileparts(tmpDir);
    
    ResultsLast2Dirs = ['..',filesep,TopDir2,filesep,TopDir1,filesep];
    
    params{IndInput}.printfile    = strrep([fName,' Overlay'],' ','_');
    
    text_annot    = cell(3,1);
    text_annot{1} = ['Structural-Image: ',structural_img];
    text_annot{2} = ['Overlay-Image: "..',ResultsLast2Dirs,fName,ext,'"'];
    text_annot{3} = ColorInfoStr; 
    
    %% make overlay & ask for saveing or not
    %% create slover object & "print" to graphics window
    if(isfield(params{IndInput},'img'))
        %% add a new window and try to leave space for the text
        params{IndInput}.figure        = spm_figure(); %new figure
        params{IndInput}.area.position = [0.005,0.005,0.99,0.95];
        params{IndInput}.area.units    = 'normalized';
        params{IndInput}.area.halign   = 'center';
        params{IndInput}.area.valign   = 'middle';
        pause(.25); %there is a strange bug that I can not nail down that disappears when waiting for a short time at this point... ??? Why??? But it works now.
        
        %% Call slover to construct SlObj for paint
        SlObj{IndInput}          = slover(params{IndInput});
        SlObj{IndInput}.printstr = [SlObj{IndInput}.printstr,' -append'];
        
        
        %% paint slices
        disp('painting...');
        paint(SlObj{IndInput});
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
    H{IndInput} = params{IndInput}.figure; %assign handle
end

%% use CheckCompEncode(CompEncodeDir) to analyse network structure from base 2 encode
CheckCompEncode(CompEncodeDir);

%% Render display?
if(isnan(RenderOverlays))
    if(strcmp('Yes',questdlg('Render all produced overlays to file?','Render Overlays?','Yes','No','Yes')))
        RenderOverlayFigures(H,NamesForRendering,[CompEncodeDir,filesep,'OverlayRendering']);
    end
else
    if(RenderOverlays)
        RenderOverlayFigures(H,NamesForRendering,[CompEncodeDir,filesep,'OverlayRendering']);
    end
end

%% done
disp(' ');
disp('Done.');
disp(' ');
end


%% subfunctions
function NIFTI_files = GetCompEncodeFiles(CompEncodeDir)
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

for IndFile = 1:length(Files)
    if(~exist(NIFTI_files{IndFile},'file'))
        error(['Can not find "',Files{IndFile},'" in "',CompEncodeDir,'".'])
    end
end

end

%% CreateColors for display
function [Colors,ColorRange,ColorInfoStr,NComp] = CreateColors(NIIfile,ColGenType,DispColors,Version)
% This function creates the colors for the overlay using distinguishable colors
% and depending on ColGenType the overlap of areas/components as white, if "SimpleEncode".

%Default setting
if(strcmp(ColGenType,'Log2Base2Encode'))
    try
        disp(['For "',ColGenType,'": Will use Version ',num2str(Version),' of color generation function. (default==1)(2==alternate)']);
    catch CATCH_Version
        disp_catch(CATCH_Version,[mfilename,'>CATCH_Version'],'CATCH_Version'); %assign & display error as warning
        Version  = 1; %FOR Log2Base2Encode: 1==make distinguishable colors as usual and interpolate inbetween + white at top&bottom; 2==make distinguishable colors as usual + white at top&bottom AND then make more distinguishable colors that are used as MidColor for mixing as before.
        disp(['For "',ColGenType,'": Will use Version ',num2str(Version),'(default) of color generation function.']);
    end
end
WaitTime = 1; %s DEFAULT time for "Color" message box staying on screen if not clicked by user

%% load NIFTI
V_Data     = spm_vol(NIIfile);
Data       = V_Data.private.dat(:);
UniqueData = unique(Data); UniqueData(UniqueData==0) = []; %zero==background is not important here (remove).

%% ColorRange
switch(ColGenType)
    case 'Base2Encode'
        disp('Skipping Base2Encode...');
        Colors       = [];
        ColorRange   = [];
        ColorInfoStr = [];
        return;
    otherwise
        ColorRange = [0.9 ceil(max(UniqueData))];
        if(strcmp(ColGenType,'Log2Base2Encode'))
            NColors = floor(max(UniqueData));
        else
            NColors = max(UniqueData);
        end
end
NComp = NColors+1; %number of "components to expect", therefore raise by one.

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
switch(ColGenType)
    case {'SimpleEncode','Log2Base2Encode','SumMap','WinnerTakeAllEncode'}
        bg = [1 0 0; 1 1 1; 0 0 0; .05 .05 .05; .1 .1 .1; .15 .15 .15; .25 .25 .25; .35 .35 .35; .45 .45 .45; .5 .5 .5; .65 .65 .65; .75 .75 .75; .85 .85 .85; .9 .9 .9]; %white black gray(s)
    otherwise
        bg = [1 1 1; 0 0 0; .05 .05 .05; .1 .1 .1; .15 .15 .15; .25 .25 .25; .35 .35 .35; .45 .45 .45; .5 .5 .5; .65 .65 .65; .75 .75 .75; .85 .85 .85; .9 .9 .9]; %white black gray(s)
end

%make colors
try
    if(NColors<=3) %yellow,blue,green as rgb
        disting_colors = [1 1 0; 0 0 1; 0 1 0]; %red is always the special case therefore no red here
        if(NColors<3)
            disting_colors((NColors+1):3,:) = [];
            if(NColors==1)
                disting_colors = [1 0 0]; %red special case because here we do not have a mixture
            end
        end
    else
        disting_colors = distinguishable_colors(NColors,bg);
    end
catch CATCH_distinguishable_colors
    disp_catch(CATCH_distinguishable_colors,[mfilename,'>CATCH_distinguishable_colors'],'CATCH_distinguishable_colors');
    disp(['NColors= ',num2str(NColors),'.']);
    return;
end
%NEW THING: SumMap gets hot colorbar but starts at blue, then green for two overlaps and then hot
if(strcmp(ColGenType,'SumMap'))
    disting_colors = hot(size(disting_colors,1));
    disting_colors(1,:) = [0 0 1]; %put green for the singular ones. (others will not be red because of bg that was set before.
    disting_colors(2,:) = [0 1 0]; %put green for two overlaps. (others will not be red because of bg that was set before.
end
    
if(Version==2)
    excluded = [disting_colors; 1 0 0; 1 1 1; 0 0 0]; %old: [disting_colors;bg];
    disting_colors2 = distinguishable_colors(NColors,excluded);
end

switch(ColGenType)
    case 'SimpleEncode'
        NTotalColors = NColors + 1; %one more for "-1" overlaps in this case
        Colors = zeros(NTotalColors,3); %init black and just change colors that are "pure areas code"
        Colors(1,:)     = [1 0 0]; %red
        Colors(2:end,:) = disting_colors;
    case 'Log2Base2Encode'
        NInterpSteps = 9; %extra interpolations between base-colors, i.e. the distinguishable colors
        NTotalColors = (NInterpSteps+1)*NColors; %make N steps between the base-colors plus the base-colors themselfes
        Colors       = ones(NTotalColors,3); %init white and just change colors that are "pure areas code", later add in the mixed colors
        for IndCol = 1:(size(disting_colors,1)+1)
            if(IndCol<=size(disting_colors,1))
                Colors(((IndCol-1)*(NInterpSteps+1))+1,:) = disting_colors(IndCol,:);
            else
                Colors(((IndCol-1)*(NInterpSteps+1))+1,:) = disting_colors(IndCol-1,:);
            end
            if(IndCol>1) %mix only starting after the first
                BottomColor = disting_colors(IndCol-1,:);
                if(IndCol<=size(disting_colors,1))
                    TopColor    = disting_colors(IndCol,:);
                else
                    TopColor    = [1 1 1]; %init this as white
                end
                if(Version==1)
                    MidColor =  (BottomColor+TopColor)/2;
                    DiffColor= TopColor-BottomColor;
                else %Version2: use a second set of colors
                    if(IndCol<=size(disting_colors,1))
                        MidColor = disting_colors2(IndCol,:);
                    else
                        MidColor = disting_colors2(IndCol-1,:);
                    end
                end
                
                %old: MixColors = [(disting_colors(IndCol-1,:)+MidColor)/2; (disting_colors(IndCol,:)+MidColor)/2];
                MixColors= zeros(NInterpSteps,3);
                for IndMix = 1:NInterpSteps
                    if(Version==2)
                        if(IndMix<=(NInterpSteps/2))
                            DiffColor= MidColor-BottomColor;
                        else
                            DiffColor= TopColor-MidColor;
                        end
                    end
                    %old: MixColors(IndMix,:) = disting_colors(IndCol-1,:)+IndMix.*DiffColor./NInterpSteps; %linear slide "upward" to the next color (in color vector space appoximated linearly, -wrong, but who cares...)
                    if(IndMix<=(NInterpSteps/2))
                        MixColors(IndMix,:) = MidColor-(ceil(NInterpSteps/2)-IndMix).*DiffColor./NInterpSteps; %linear slide "upward" to the next color (in color vector space appoximated linearly, -wrong, but who cares...)
                    else
                        MixColors(IndMix,:) = MidColor+IndMix.*DiffColor./(2.1*NInterpSteps); %linear slide "upward" to the next color (in color vector space appoximated linearly, -wrong, but who cares...)
                    end
                    %check NB: these "check" are a bit bad right now, maybe we can fix them in a future version.
                    if(any(MixColors(IndMix,:)>1))
                        MixColors(IndMix,MixColors(IndMix,:)>1) = 1.*ones(size(MixColors(IndMix,MixColors(IndMix,:)>1)))-(MixColors(IndMix,MixColors(IndMix,:)>1)-1.*ones(size(MixColors(IndMix,MixColors(IndMix,:)>1))));
                        disp([num2str(IndCol),',',num2str(IndMix),': >1']);
                    end
                    if(any(MixColors(IndMix,:)<0))
                        MixColors(IndMix,MixColors(IndMix,:)<0) = -MixColors(IndMix,MixColors(IndMix,:)<0);
                        disp([num2str(IndCol),',',num2str(IndMix),': <0']);
                        if(any(MixColors(IndMix,:)>1))
                            MixColors(IndMix,MixColors(IndMix,:)>1) = MixColors(IndMix,MixColors(IndMix,:)>1)-1.*ones(size(MixColors(IndMix,MixColors(IndMix,:)>1)));
                            disp('...and now: >1');
                        end
                    end
                end
                if(IndCol<=size(disting_colors,1))
                    Colors((((IndCol-2)*(NInterpSteps+1))+2):((IndCol-1)*(NInterpSteps+1)),:)                     = MixColors;
                else
                    Colors((((IndCol-2)*(NInterpSteps+1))+2):(((IndCol-2)*(NInterpSteps+1))+2)+size(MixColors,1)-1,:) = MixColors;
                end
            end
        end
        if(size(Colors,1)>NTotalColors)
            disp('Ups, we have raised the roof! Lowering it again.');
            Colors = Colors(1:NTotalColors,:); %remove extra
        end
        Colors = circshift(Colors,max([1,floor(NInterpSteps/1.85)])); %shift to arrange colors more neatly
    otherwise
        NTotalColors = NColors;
        Colors = disting_colors;
end

%% Display colors
if(DispColors)
    AllColorsExample  = repmat([size(Colors,1):-1:1]',1,10); %example for plotting
    switch(ColGenType)
        case {'SimpleEncode','Base2Encode','UniqueOnly'}
            figure; imagesc(AllColorsExample);  colormap(Colors);         title(['All possible colors (',num2str(NTotalColors),') including mixing for area overlaps (ColGenType: "',ColGenType,'")']); colorbar; axis('off');
        case 'Log2Base2Encode'
            BaseColorsExample  = repmat([size(disting_colors,1):-1:1]',1,10); %example for plotting
            figure; imagesc(AllColorsExample);   colormap(Colors);         title(['All possible colors (',num2str(NTotalColors),') including mixing for area overlaps (ColGenType: "',ColGenType,'") for ',num2str(NColors),' Components.']); colorbar; axis('off');
            figure; imagesc(BaseColorsExample);  colormap(disting_colors); title(['All possible BASIC colors (',num2str(NColors),') WITHOUT mixing for area overlaps (ColGenType: "',ColGenType,'")']); colorbar; axis('off');                            
        otherwise
            figure; imagesc(AllColorsExample);  colormap(Colors);         title(['All possible colors (',num2str(NTotalColors),') [ColGenType: "',ColGenType,'"]']); colorbar; axis('off');
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
end

%% info string
if(NTotalColors>1)
    ColorInfoStr = [num2str(NTotalColors),' Colors generated for ColGenType: "',ColGenType,'"'];
else
    ColorInfoStr = ['Color generated for ColGenType: "',ColGenType,'"'];
end
    

end