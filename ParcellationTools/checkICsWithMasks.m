function [ResultsStruct,H,AX] = checkICsWithMasks(ThreshNIIpaths,MasksPaths,SaveNormICs)
% This function can check the thresholded(!) ICs, i.e. NIFTIs-files where "0" indicates not significant  
% and other values t- or z-scores, NOT P-values (sorry for that) against input images (assumed to be masks).
% The function outputs a listing of ICs associated with the masks and returns a structure containing the listing and further information.
%
% The listing is based on a kind of probability measure that is determined in three steps.
%   1.Take the absolute of the statistic values in the thresholded NIFTI.
%   2.Sum all those values and then divide each value by this sum,
%     i.e. normalize the sum to be equal to "1" in total.
%   3.Determine the sum of the normalized values for each mask,
%     i.e. determine how much each mask contributes to the total sum.
%
% NB: This will give a kind of probability value that captures overlap and
%     significance amplitude for each mask VS the total brain, and should
%     therefore allow a judgement of how much a particular mask does contribute, 
%     in terms of the size of overlap and the amplitude of the significance value. 
%
%
%Inputs:
%        ThreshNIIpaths & MasksPaths have to be cellstrings (or convertible char-arrays) pointing to the paths to the ICs and Masks, respectively.
%Outputs:
%        ResultsStruct contains all the information and Results of checkICsWithMasks processing using inputs ThreshNIIpaths and MasksPaths
%        H & AX are the figure handles and axis handles, respectively, that were used during plotting.
%
%
%
%Usage:
%       [ResultsStruct,H,AX] = checkICsWithMasks(ThreshNIIpaths,MasksPaths); 
%       [ResultsStruct,H,AX] = checkICsWithMasks(); %select input data manually
%        ResultsStruct       = checkICsWithMasks(); %NO PLOT (only one output==ResultsStruct) & select input data manually
%       [ResultsStruct,H,AX] = checkICsWithMasks(ResultsStruct); %just use the plotting functionality with an existing ResultsStruct.
%                              checkICsWithMasks(ResultsStruct); %IN THIS CASE PLOTTING IS ALWAYS DONE EVEN IF NO OUTPUT IS SPECIFIED. Just use the plotting functionality with an existing ResultsStruct.
% 
%
%V1.6
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.6(26.10.2016): If empty selection with spm_select is made for MasksPaths, then it will be replaced with ThreshNIIpaths.  V1.5(11.10.2016): Add Dice Coefficient and SQRT of the Product of both as well. V1.0: (08.10.2016): initial implementation based on previous tests

%% Check inputs
%SPECIAL CASE, INPUT IS A STRUCT
%--> assume it is the results struct, then just display its contents and return
if(nargin==1)
    if(isstruct(ThreshNIIpaths))
        ResultsStruct = ThreshNIIpaths; %assume it is the results structure
        [H,AX] = DisplayResults(ResultsStruct);
        disp('Done with plotting ResultsStruct');
        disp(' ');
        return;
    else
        if(ischar(ThreshNIIpaths)&&length(cellstr(ThreshNIIpaths))==1) %check if it is a path to a mat-file
            [BaseDir,FName,Ext] = fileparts(ThreshNIIpaths);
            if(strcmp(Ext,'.mat'))
                disp(['First input seems to be a path to a *.mat file ("',FName,Ext,'" from directory "',BaseDir,'")']);
                disp( 'Will try to find the ResultsStruct and then load it for plotting it.');
                data = whos(matfile(ThreshNIIpaths));
                if(any(strcmp({data.name},'ResultsStruct')))
                    load(ThreshNIIpaths);
                    [H,AX] = DisplayResults(ResultsStruct);
                    disp('Done with plotting ResultsStruct');
                    disp(' ');
                    return;
                else
                    error(['Could not find "ResultsStruct" in file "',FName,Ext,'" from directory "',BaseDir,'"! Check input.']);
                end
            end
            clear BaseDir FName Ext
        end                
    end
end

%ThreshNIIpaths
if(~exist('ThreshNIIpaths','var'))
    ThreshNIIpaths = cellstr(spm_select(Inf,'image','Select Thresholded ICs...',{},pwd,'^thresh_zstat'));
else
    if(isempty(ThreshNIIpaths))
        ThreshNIIpaths = cellstr(spm_select(Inf,'image','Select Thresholded ICs...',{},pwd,'^thresh_zstat'));
    else
        if(~iscellstr(ThreshNIIpaths))
            if(~ischar(ThreshNIIpaths))
                error('ThreshNIIpaths must be a cellstring! (or convertible!)');
            else
                ThreshNIIpaths = cellstr(ThreshNIIpaths);
            end
        end
    end
end

%MasksPaths
if(~exist('MasksPaths','var'))
    MasksPaths = cellstr(spm_select(Inf,'image','Select Masks...'));
else
    if(isempty(MasksPaths))
        MasksPaths = cellstr(spm_select(Inf,'image','Select Masks...'));
    else
        if(~iscellstr(MasksPaths))
            if(~ischar(MasksPaths))
                error('MasksPaths must be a cellstring! (or convertible!)');
            else
                MasksPaths = cellstr(MasksPaths);
            end
        end
    end
end
if(length(MasksPaths)==1&&isempty(MasksPaths{1})) %assume user wants to compare ThreshNIIpaths with themselves
    MasksPaths = ThreshNIIpaths;
end

%SaveNormICs
if(~exist('SaveNormICs','var'))
    SaveNormICs = 0;
else
    if(isempty(SaveNormICs))
        SaveNormICs = 0; %just to be save.
    end
end

%% Do a quick check of the masks and ICs to see if the dimensions are the same.
%% Also sort ICs!!! (if possible)
[ok,ICnums] = findICnums(ThreshNIIpaths);
if(~ok)
    ICnums = 1:length(ThreshNIIpaths);
else
    disp('Sorting input ICs...');
    [ICnums,sortInds] = sort(ICnums);
    ThreshNIIpaths    = ThreshNIIpaths(sortInds); %InputFiles must be a cellstr otherwise this doesn't work therefore the hassle above when checking the inputs.
    clear sortInds
end
VICs  = spm_vol(ThreshNIIpaths);
VMasks= spm_vol(MasksPaths);
if(~spm_check_orientations([VICs{:},VMasks{:}]))
    if(~spm_check_orientations([VICs{:}]) && ~spm_check_orientations([VMasks{:}]))
        error('Dimensions of ICs, as well as, Masks are not all the same! (Differences between ICs detected, and also differences between Masks detected!)');
    else
        if(~spm_check_orientations([VICs{:}]))
            error('Dimensions of ICs are not the all same!');
        elseif(~spm_check_orientations([VMasks{:}]))
            error('Dimensions of Masks are not the all same!');
        else
            error('Dimensions of ICs and Masks are not the same! (However, ICs are internally consistent, and separately Masks seem to be internally consistent.)');
        end
    end
end

%% setup important constants and templates for further use
Factor  = 100; %normalize sums to this Factor 
NICs    = length(ThreshNIIpaths);
NMasks  = length(MasksPaths);
VTmpIC  = VICs{1}; %take the first one as an example
VTmpMask= VMasks{1};     %take the first one as an example

%% load ThreshNII & MaskNII data
ICdata4D = zeros(NICs,    VTmpIC.dim(1),  VTmpIC.dim(2),  VTmpIC.dim(3));
Masks4D  = zeros(NMasks,VTmpMask.dim(1),VTmpMask.dim(2),VTmpMask.dim(3));
for IndIC = 1:NICs
    ICdata4D(IndIC,:,:,:)  = VICs{IndIC}.private.dat(:,:,:,VICs{IndIC}.n(1));
end
for IndMask = 1:NMasks
    Masks4D(IndMask,:,:,:) = VMasks{IndMask}.private.dat(:,:,:,VMasks{IndMask}.n(1))~=0;
end

%% get total mask based on non-zero values of all ICs
ICsMask3D  = squeeze(sum(double(ICdata4D~=0),1));
Masks3D    = squeeze(sum(double( Masks4D~=0),1));
if(any(Masks3D(:)>1))
    MasksHaveOverlap = 1;
    if(nargout>1)
        Hwait = helpdlg({'WARNING:';...
            'Some of the masks are overlapping!';...
            'This will offset the sums over the all maps,';
            'such that it will be larger than 1!';
            ' ';
            'Be careful with the interpretation!!!'},'WARNING!');
        uiwait(Hwait);
    else
        disp(' ');
        disp('WARNING: Some of the masks are overlapping!');...
        disp('         This will offset the sums over the all maps, such that it will be larger than 1!');
        disp(' ');
        disp('         Be careful with the interpretation!!!');
        disp(' ');
    end
else
    MasksHaveOverlap = 0;
end
%The TOTAL MASK!
TotalMask3D   = ICsMask3D|Masks3D; 
NVoxTotalMask = length(find(TotalMask3D(:)~=0));

%% get IC data and all masks in the total mask
ICdata2D = zeros(NICs,NVoxTotalMask); %The ICs in the total mask
for IndIC = 1:NICs
    CurrData3D = squeeze(ICdata4D(IndIC,:,:,:));
    ICdata2D(IndIC,:) = CurrData3D(TotalMask3D(:)~=0);
end

Masks2D  = zeros(NMasks,NVoxTotalMask); %Each of the Masks in the total mask
NVoxMasks= zeros(NMasks,1);
for IndMask = 1:NMasks
    CurrData3D = Masks4D(IndMask,:,:,:);
    NVoxMasks(IndMask) = length(find(CurrData3D(:)~=0));
    Masks2D(IndMask,:) = CurrData3D(TotalMask3D(:)~=0);
end

%% normalize absolute values of each IC
absICdata2D = abs(ICdata2D);
NormVals    = sum(absICdata2D,2);
normICs     = Factor.*absICdata2D./repmat(NormVals,1,NVoxTotalMask);

%% SaveNormICs?
if(SaveNormICs)
    disp('Will save the sum-normalized IC maps with prefix "SumNorm_"...');
    normICdata4D = Factor.*abs(ICdata4D)./repmat(NormVals,1,size(ICdata4D,2),size(ICdata4D,3),size(ICdata4D,4));
    for IndIC = 1:NICs
        CurrData = squeeze(normICdata4D(IndIC,:,:,:));
        CurrData = CurrData./repmat(max(CurrData(:)),size(CurrData));
        normICdata4D(IndIC,:,:,:) = CurrData;
        
        WriteNormIC(VICs{IndIC},ThreshNIIpaths{IndIC},'SumNorm_',CurrData);
    end
    
    %save also a copy of the sum over all norm ICs AND mean over all norm ICs
    SumAllNormICs  = squeeze(sum( normICdata4D,1));
    MeanAllNormICs = squeeze(mean(normICdata4D,1));
    
    Vtemplate = rmfield(VICs{1},'private');
    if(Vtemplate.n(1)~=1)
        Vtemplate.n(1) = 1; %this can not be a 4D file
    end
    OutDirSumAllNormICs = spm_select(1,'dir','Select an output directory for the the sum over all & mean over all sum-normlized thresholded ICs...');
    if(~isempty(OutDirSumAllNormICs))
        %SumAllNormICs
        Vtemplate.fname = [OutDirSumAllNormICs,filesep,'SumAllNormICs.nii'];
        WriteNormIC(Vtemplate,Vtemplate.fname,'',SumAllNormICs);
        
        %MeanAllNormICs
        Vtemplate.fname = [OutDirSumAllNormICs,filesep,'MeanAllNormICs.nii'];
        WriteNormIC(Vtemplate,Vtemplate.fname,'',MeanAllNormICs);
    else
        disp('No output directory selected for saving SumAllNormICs.nii & MeanAllNormICs.nii, will skip saving.');
    end
    
    
end

%% free some memory
if(SaveNormICs)
    clear ICdata4D Masks4D normICdata4D SumAllNormICs MeanAllNormICs
else
    clear ICdata4D Masks4D
end

%% MaskNames & also get size of Mask in the overall Mask --> use this for exclusion later.
[BaseDirs,MaskNames,Ext] = cellfun(@fileparts,MasksPaths,'UniformOutput',false);
MaskNames = regexprep(MaskNames,'_',' ');
clear BaseDirs Ext

%% for each normIC check the sum of each mask --> make a matrix "SumContributions" IC#-x-NMasks
SumContributions = zeros(NICs,NMasks);
for IndMask = 1:NMasks
    SumContributions(:,IndMask) = sum(normICs(:,Masks2D(IndMask,:)~=0),2);
end

%% also for each normIC after creating SumContributions Matrix make a MAX-SCALED verions of this, 
%  i.e. --> make matrix "MaxScaleSumContributions" where the matrix SumContributions was divided by the maximum of each row (i.e. best/highest value mask)
MaxScaleSumContributions = SumContributions./repmat(max(SumContributions,[],2),1,NMasks);

%% calculate Dice Coefficients
% DiceCoeff = zeros(NICs,NMasks);
% for IndMask = 1:NMasks
%     for IndIC = 1:NICs
%         DiceCoeff(IndIC,IndMask) = 2.*sum(double((normICs(IndIC,:)~=0) & (Masks2D(IndMask,:)~=0)))./(sum(double(normICs(IndIC,:)~=0))+sum(double(Masks2D(IndMask,:)~=0)));
%     end
% end    
A = double(normICs~=0);
B = double(Masks2D~=0);
OverlapSize= A*B';
DiceCoeff  = 2.*(A*B')./(repmat(sum(A,2),1,size(B,1))+repmat(sum(B,2)',size(A,1),1));

%% sqrt(Product)
ProductCoeff = sqrt(MaxScaleSumContributions.*DiceCoeff);

%% create Struct-Array ICassignments with three fields SortedVals & SortInds & Name which is either SumContributions or DiceCoeff
%  These two things can be plotted in a bar graph with extra text to indicate the MaskNames FOR EACH IC, selected by the user.
%  Before choosing the IC the user should also get a list of ICs that only have zero values for all the Masks --> i.e. not overlapping the masks at all!
ICassignments(1).Name       = 'SumContributions';
ICassignments(1).SortedVals = zeros(NICs,NMasks);
ICassignments(1).SortInds   = zeros(NICs,NMasks);

ICassignments(2).Name       = 'DiceCoeff';
ICassignments(2).SortedVals = zeros(NICs,NMasks);
ICassignments(2).SortInds   = zeros(NICs,NMasks);

ICassignments(3).Name       = 'ProductCoeff';
ICassignments(3).SortedVals = zeros(NICs,NMasks);
ICassignments(3).SortInds   = zeros(NICs,NMasks);

ICassignments(4).Name       = 'OverlapSize';
ICassignments(4).SortedVals = zeros(NICs,NMasks);
ICassignments(4).SortInds   = zeros(NICs,NMasks);
for IndAssign = 1:length(ICassignments)
    if(IndAssign==1)
        Data = SumContributions;
    elseif(IndAssign==2)
        Data = DiceCoeff;
    elseif(IndAssign==3)
        Data = ProductCoeff;
    else
        Data = OverlapSize;
    end
    for IndIC = 1:NICs
        [SortedVals,SortInds]  = sort(Data(IndIC,:),'descend');
        if(any(SortedVals(:)~=0))
            ICassignments(IndAssign).SortedVals(IndIC,:) = SortedVals;
            ICassignments(IndAssign).SortInds(IndIC,:)   = SortInds;
        end
    end
end

%% Pack data into ResultsStruct such that it can be reused.
ResultsStruct = struct('ThreshNIIpaths',{ThreshNIIpaths},'MasksPaths',{MasksPaths},'ICnums',ICnums,'NICs',NICs,'NMasks',NMasks,'MaskNames',{MaskNames},'NVoxMasks',NVoxMasks,'MasksHaveOverlap',MasksHaveOverlap,'normICs',normICs,'NormVals',NormVals,'Factor',Factor,'SumContributions',SumContributions,'MaxScaleSumContributions',MaxScaleSumContributions,'DiceCoeff',DiceCoeff,'ProductCoeff',ProductCoeff,'OverlapSize',OverlapSize,'ICassignments',ICassignments);

%% Do Display if outputs H & AX are set.
if(nargout>1)
    [H,AX] = DisplayResults(ResultsStruct);
end

%% Done.
disp('DONE.');
disp(' ');

end

%% Subfunctions
%% DisplayResults
function [H,AX] = DisplayResults(ResultsStruct)
% This function displays the Data in ResultsStruct

%% display SumContributions matix and MaxScaleSumContributions side by side
%  & add names of masks to display (x-axis of first matrix get even# MaskNames and x-axis of second Matirx gets odd# MaskNames to keep it readable, hopefully)
LimStr = ['[0,',num2str(ResultsStruct.Factor),']'];
ChoiceLimits = questdlg({'Which limits shall be used for the matrix plot of SumContributions?';['Minimum to maximum, OR fixed ',LimStr,', OR manually set?']},'Limits',LimStr,'Min2Max','Manual','Min2Max');
switch(ChoiceLimits)
    case LimStr
        Limits = [0,ResultsStruct.Factor];
    case 'Min2Max'
        Limits = [min(ResultsStruct.SumContributions(:)),max(ResultsStruct.SumContributions(:))];
    otherwise
        AnswerLimits = inputdlg({'Limits for SumContributions?'},'Limits',1,{LimStr});
        Limits = str2num(AnswerLimits{1});
end

MainFig = ResultsStruct.NICs+10;
while(ishandle(MainFig)) %find a figure above NICs 
    MainFig = MainFig + 1;
end
H{1} = figure(MainFig); clf;
% AX{1} = subplot(1,4,1); 
imagesc(ResultsStruct.SumContributions,Limits); title({'SumContributions IC#-x-NMasks';['Contributions to total sum==',num2str(ResultsStruct.Factor)]});
set(AX{1},'XTick',2:2:ResultsStruct.NMasks); %set which ticks to set
set(AX{1},'XTickLabel',ResultsStruct.MaskNames(2:2:end)) % set tick labels; here even numbered names
set(AX{1},'XTickLabelRotation',45);          %rotate labels by 45°

H{1} = figure(MainFig+1); clf;
% AX{2} = subplot(1,4,2); 
imagesc(ResultsStruct.MaxScaleSumContributions,[0 1]); title({'MaxScaleSumContributions IC#-x-NMasks';'MaxScaled: highlighting contributions.'});
set(AX{2},'XTick',1:2:ResultsStruct.NMasks); %set which ticks to set
set(AX{2},'XTickLabel',ResultsStruct.MaskNames(1:2:end)) % set tick labels; here odd numbered names
set(AX{2},'XTickLabelRotation',45);          %rotate labels by 45°

H{1} = figure(MainFig+2); clf;
% AX{3} = subplot(1,4,3);
imagesc(ResultsStruct.ProductCoeff); title('Product Coefficient');
set(AX{3},'XTick',2:2:ResultsStruct.NMasks);      %set which ticks to set
set(AX{3},'XTickLabel',ResultsStruct.MaskNames(2:2:end)) %set tick labels; here odd numbered names
set(AX{3},'XTickLabelRotation',45);             %rotate labels by 45°

H{1} = figure(MainFig+3); clf;
% AX{4} = subplot(1,4,4);
imagesc(ResultsStruct.DiceCoeff); title('Dice Coefficient');
set(AX{4},'XTick',1:2:ResultsStruct.NMasks);      %set which ticks to set
set(AX{4},'XTickLabel',ResultsStruct.MaskNames(1:2:end)) %set tick labels; here odd numbered names
set(AX{4},'XTickLabelRotation',45);             %rotate labels by 45°

H{1} = figure(MainFig+4); clf;
% AX{4} = subplot(1,4,4);
imagesc(ResultsStruct.OverlapSize); title('Overlap Size (i.e. cluster size and overlaps)');
set(AX{4},'XTick',1:2:ResultsStruct.NMasks);      %set which ticks to set
set(AX{4},'XTickLabel',ResultsStruct.MaskNames(1:2:end)) %set tick labels; here odd numbered names
set(AX{4},'XTickLabelRotation',45);             %rotate labels by 45°

%% THRESHOLDED PLOT
H{2} = figure(MainFig+5); clf;
AX{4+1} = subplot(1,4,1); imagesc(ResultsStruct.SumContributions.*double(ResultsStruct.SumContributions>(ResultsStruct.Factor.*0.1)),Limits); title({'SumContributions IC#-x-NMasks > 10%';['Contributions to total sum==',num2str(ResultsStruct.Factor)]});
set(AX{4+1},'XTick',2:2:ResultsStruct.NMasks); %set which ticks to set
set(AX{4+1},'XTickLabel',ResultsStruct.MaskNames(2:2:end)) % set tick labels; here even numbered names
set(AX{4+1},'XTickLabelRotation',45);          %rotate labels by 45°

AX{4+2} = subplot(1,4,2); imagesc(ResultsStruct.MaxScaleSumContributions.*double(ResultsStruct.MaxScaleSumContributions>0.1),[0 1]); title({'MaxScaleSumContributions IC#-x-NMasks > 10%';'MaxScaled: highlighting contributions.'});
set(AX{4+2},'XTick',1:2:ResultsStruct.NMasks); %set which ticks to set
set(AX{4+2},'XTickLabel',ResultsStruct.MaskNames(1:2:end)) % set tick labels; here odd numbered names
set(AX{4+2},'XTickLabelRotation',45);          %rotate labels by 45°

AX{4+3} = subplot(1,4,3); imagesc(ResultsStruct.ProductCoeff.*double(ResultsStruct.ProductCoeff>0.1)); title('Product Coefficient > 10%');
set(AX{4+3},'XTick',2:2:ResultsStruct.NMasks);      %set which ticks to set
set(AX{4+3},'XTickLabel',ResultsStruct.MaskNames(2:2:end)) %set tick labels; here odd numbered names
set(AX{4+3},'XTickLabelRotation',45);             %rotate labels by 45°

AX{4+4} = subplot(1,4,4); imagesc(ResultsStruct.DiceCoeff.*double(ResultsStruct.DiceCoeff>0.1)); title('Dice Coefficient > 10%');
set(AX{4+4},'XTick',1:2:ResultsStruct.NMasks);      %set which ticks to set
set(AX{4+4},'XTickLabel',ResultsStruct.MaskNames(1:2:end)) %set tick labels; here odd numbered names
set(AX{4+4},'XTickLabelRotation',45);             %rotate labels by 45°


%% display bar graph for a USER CHOICE IC with contribution and Mask Names indicated.
%  Before choosing the IC the user should also get a list of ICs that only have zero values for all the Masks --> i.e. not overlapping the masks at all!
ICsExamined  = []; %init empty
ICsRemaining = (1:ResultsStruct.NICs)'; %init
KeepPlotting = 1;
while(KeepPlotting)
    %select an IC
    [SelInds,OK] = listdlg('ListString',cellstr(num2str(ResultsStruct.ICnums(ICsRemaining))),'SelectionMode','multiple','Name','ICselect','PromptString','Select an IC for plotting:','CancelString','Quit');
    if(~OK)
        disp('Quit by user choice.');
        break;
    end
    SelIC = ICsRemaining(SelInds); %get transformation right.
    
    ICsExamined = unique([ICsExamined; SelIC(:)]); %add to list
    for Ind = 1:length(SelIC)
        if(isempty(ICsRemaining))
            break;
        end
        ICsRemaining(ICsRemaining==SelIC(Ind)) = []; %remove this one from list.
        H{end+1} = figure(); clf; %get figure
        for IndAssign = 1:length(ResultsStruct.ICassignments)
            %plot bargraph with names added.
            if(IndAssign==1)
                CurrFactor= ResultsStruct.Factor;
            else
                CurrFactor= 1;
            end
            CurrName  = ResultsStruct.ICassignments(IndAssign).Name;
            CurrVals  = ResultsStruct.ICassignments(IndAssign).SortedVals(SelIC(Ind),:);
            CurrInds  = ResultsStruct.ICassignments(IndAssign).SortInds(SelIC(Ind),:);
            ReleventMasksIC = ResultsStruct.MaskNames(CurrInds);
            
            if(all(CurrVals==0))
                Hwait = helpdlg({['IC ',num2str(ResultsStruct.ICnums(SelIC(Ind))),' has NO overlap with any of the Masks!']; 'I.e. all SumContributions are zero in all Masks.'; ' '; 'That is interesting too, is it not?'},['No Overlap IC ',num2str(SelIC(Ind))]);
                uiwait(Hwait);
            else
                if(ResultsStruct.MasksHaveOverlap) %in this case CumSum doesn't make sense / is misleading.
                    IndCumSum80 = []; %DON'T USE! %index where we have 80% cumulative sum (or empty if not exist)
                    IndCumSum50 = []; %DON'T USE! %index where we have 50% cumulative sum (or empty if not exist)
                else
                    CumSumIC = cumsum(CurrVals);
                    IndCumSum80 = find(CumSumIC>=(CurrFactor.*0.80),1,'first'); %index where we have 80% cumulative sum (or empty if not exist)
                    IndCumSum50 = find(CumSumIC>=(CurrFactor.*0.50),1,'first'); %index where we have 50% cumulative sum (or empty if not exist)
                end
                
                %plot current data into subplot
                subplot(1,length(ResultsStruct.ICassignments),IndAssign);
                bar(1:length(CurrVals),CurrVals); hold on
                plot(1:length(CurrVals),CurrFactor.*0.05.*ones(1,length(CurrVals)),'k--');
                plot(1:length(CurrVals),CurrFactor.*0.10.*ones(1,length(CurrVals)),'k-');
                if(~isempty(IndCumSum80))
                    plot(repmat(IndCumSum80+0.5,1,length(0:0.1:1)),0:0.1:1,'r--');
                    legend({'Contributions in each Mask';'5% contribution';'10% contribution';'80% of IC covered'});
                else
                    if(~isempty(IndCumSum50))
                        plot(repmat(IndCumSum50+0.5,1,length(0:0.1:1)),0:0.1:1,'r--');
                        legend({'Contributions in each Mask';'5% contribution';'10% contribution';'50% of IC covered'});
                    else
                        legend({'Contributions in each Mask';'5% contribution';'10% contribution'});
                    end
                end
                set(gca,'XTick',1:length(CurrVals));   %which ticks to set
                set(gca,'XTickLabel',ReleventMasksIC); %set tick labels;
                set(gca,'XTickLabelRotation',45);      %rotate labels by 45°
                title([CurrName,' for IC ',num2str(ResultsStruct.ICnums(SelIC(Ind))),' in decending order.']);
            end
        end
    end
    if(~isempty(ICsRemaining) && (length(ICsExamined)~=ResultsStruct.NICs))
        %Plot another one?
        if(strcmp('No',questdlg('Look at another IC?','Another?','Yes','No','Yes')))
            disp('User chose to quit plotting.');
            KeepPlotting = 0;
        end
    else
        disp('All ICs have been plotted.');
        KeepPlotting = 0;
    end
end
disp(' ');

end


%% WriteNormIC(VICs{IndIC},ThreshNIIpaths{IndIC},'SumNorm_',squeeze(normICdata4D(IndIC,:,:,:)))
function [Vout,OutputFName] = WriteNormIC(Vtmp,OrgNIIpath,Prefix,ICdata3D)
% This function simply writes out the data to a path and adds a prefix to the filename

%% create new filepath
[BaseDir,FName,Ext] = fileparts(OrgNIIpath);
if(strcmpi(Ext,'.nii'))
    OutputFName= [Prefix,FName,Ext];
else
    OutputFName= [Prefix,FName,'.nii'];
end
OutputPath = [BaseDir,filesep,OutputFName];

%% treat special case in which this function is used to write new data and BaseDir does not exist
if(~exist(BaseDir,'dir'))
    disp(['Directory "',BaseDir,'" does not exist, will create it...']);
    mkdir(BaseDir);
end

%% Create Output Volume Struct
if(isfield(Vtmp,'private')) %just to be sure.
    Vout = rmfield(Vtmp,'private'); %just to be sure.
else
    Vout = Vtmp;
end
Vout.fname = OutputPath;
if(Vout.dt(1)<16)
    Vout.dt(1) = 16;
end

%% Write out NIFTI
disp(['Writing "',OutputFName,'" to directory "',BaseDir,'"...']);
Vout = spm_write_vol(Vout,ICdata3D);


end
