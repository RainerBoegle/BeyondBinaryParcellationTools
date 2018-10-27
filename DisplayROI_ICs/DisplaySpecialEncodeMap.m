function [H,SlObj,params,SliceNumbers,SpecialEncodePath,InfoStruct,ColorsUniques,ColorsOverlaps,AllColors] = DisplaySpecialEncodeMap(SpecialEncodePath,ColorsOverlaps,SliceNumbers)
% This function can be used for displaying the SpecialEncodeMap,
% or creating and displaying it.
%
%
%Usage:
%       [H,SlObj,params,SliceIndices,SpecialEncodePath,InfoStruct,ColorsUniques,ColorsOverlaps,AllColors] = DisplaySpecialEncodeMap(SpecialEncodePath,ColorsOverlaps,SliceNumbers);
%       [H,SlObj,params,SliceIndices,SpecialEncodePath,InfoStruct,ColorsUniques,ColorsOverlaps,AllColors] = DisplaySpecialEncodeMap();
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (06.10.2016): initial implementation based on previous tests

%% Check inputs
%SpecialEncodePath
if(~exist('SpecialEncodePath','var'))
    [SpecialEncodePath,InfoStruct,CombineOverlapNumbers] = SelectOrCreate([]);
elseif(isempty(SpecialEncodePath))
    [SpecialEncodePath,InfoStruct,CombineOverlapNumbers] = SelectOrCreate([]);
else
    if(iscellstr(SpecialEncodePath))
        SpecialEncodePath = SpecialEncodePath{1};
        if(~ischar(SpecialEncodePath))
            error('SpecialEncodePath is not a string indicating the path to the SpecialEncodeMap.');
        end
    else
        if(~ischar(SpecialEncodePath))
            error('SpecialEncodePath is not a string indicating the path to the SpecialEncodeMap.');
        end
    end
    [SpecialEncodePath,InfoStruct,CombineOverlapNumbers] = SelectOrCreate(SpecialEncodePath);
end

%ColorsOverlaps
TemplateColors = [1,0,0]; %red %old: [0,1,0]; %green %Old: [0,0,1]; %blue OLD: 
if(~exist('ColorsOverlaps','var'))
    ColorsOverlaps = CreateColorsOverlaps(CombineOverlapNumbers,TemplateColors);
elseif(isempty(ColorsOverlaps))
    ColorsOverlaps = CreateColorsOverlaps(CombineOverlapNumbers,TemplateColors);
else
    if(size(ColorsOverlaps,2)~=3)
        error('ColorsOverlaps must be of size N-x-3 and contain RGB-values for indicating colors for the Overlaps.');
    else
        if(any(ColorsOverlaps(:)<0)||any(ColorsOverlaps(:)>1))
            error('ColorsOverlaps RGB-values must from in the interval [0, 1]!');
        end
    end
end

%SliceNumbers
if(~exist('SliceNumbers','var'))
    SliceNumbers = [];
elseif(~isempty(SliceNumbers))
    disp('Using input slice numbers in mm coordinates...');
else
    disp('Will determine Slice numbers automatically...');
end

%% create all colors
%just making sure that the gray and black tones are ignored by distinguishable_colors! ;)
bg = [0 0 0; .05 .05 .05; .1 .1 .1; .15 .15 .15; .25 .25 .25; .35 .35 .35; .45 .45 .45; .5 .5 .5]; %white black gray(s)
ColorsUniques = distinguishable_colors(size(InfoStruct.ICcorrespondence,1),[ColorsOverlaps;bg]); %make colors for the components and exclude the colors for the Overlaps
AllColors     = [ColorsOverlaps;ColorsUniques];

%% Display with DisplayClusters
if(isempty(SliceNumbers))
    [H,SlObj,params,SliceNumbers,AllColors] = DisplayClusters('Cluster',SpecialEncodePath,[],AllColors);
else
    [H,SlObj,params,SliceNumbers,AllColors] = DisplayClusters(SliceNumbers,SpecialEncodePath,[],AllColors);
end

%% Done.
disp('DONE.');
disp(' ');

end

%% subfunction
%% Select Or Create SpecialEncodeMap?
function [SpecialEncodePath,InfoStruct,CombineOverlapNumbers] = SelectOrCreate(SpecialEncodePath)
% Select or create SpecialEncodeMap?
if(~isempty(SpecialEncodePath))
    CompEncodeDir = fileparts(SpecialEncodePath);
    load(spm_select('FPList',CompEncodeDir,'^SpecialEncodeInfoStruct.mat'));
    CombineOverlapNumbers = InfoStruct.CombineOverlapNumbers;
else
    %ask user.
    if(strcmp('Select',questdlg('Select SpecialEncodeMap or create it?','Select?','Select','Create','Select')))
        CompEncodeDir     = spm_select(1,'dir','Select Component Encode Directory...');
        SpecialEncodePath = spm_select('FPList',CompEncodeDir,'^SpecialEncodeMap.nii');
        load(spm_select('FPList',CompEncodeDir,'^SpecialEncodeInfoStruct.mat'));
        CombineOverlapNumbers = InfoStruct.CombineOverlapNumbers;
    else
        %Create
        CompEncodeDir = spm_select(1,'dir','Select Component Encode Directory...');
        [Hfig,ExpBinCenters] = PlotSumMapHist(spm_select('FPList',CompEncodeDir,'^SumMap.nii'));
        H = helpdlg({'1.Examine the SumMap and the histogram of the SumMap.'; '2.indicate the overlap numbers that should be combined.'},'How to...');
        uiwait(H);
        answerNCombinations = inputdlg({'How many combinations of overlaps?'},'Number of overlaps to be combined?',1,{'2'});
        NCombinations = str2num(answerNCombinations{1});
        CombineOverlapNumbers = cell(NCombinations,1);
        OverlapsRemaining = ExpBinCenters; %init
        OverlapsRemaining(OverlapsRemaining==1) = [];
        for Ind = 1:NCombinations
            answerCurrCombOverlap = inputdlg({['Overlap Combination ',num2str(Ind),':']},'Combinations?',1,{['[',JoinStr(OverlapsRemaining,','),']']});
            CurrCombOverlap = str2num(answerCurrCombOverlap{1});
            CombineOverlapNumbers{Ind} = CurrCombOverlap;
            for IndRem = 1:length(CurrCombOverlap)
                OverlapsRemaining(OverlapsRemaining==CurrCombOverlap(IndRem)) = [];
            end
            if(isempty(OverlapsRemaining))
                if(Ind~=NCombinations)
                    disp('Ending the overlap combinations, all possible values were set already.');
                    CombineOverlapNumbers(Ind+1:end) = [];
                    break;
                end
            end
        end
        [SpecialEncodePath,InfoStruct] = CreateSpecialEncodeMap(CompEncodeDir,CombineOverlapNumbers);
    end
end

end

%% CreateColorsOverlaps
function ColorsOverlaps = CreateColorsOverlaps(CombineOverlapNumbers,TemplateColor)
%ask user for the colors for the overlaps
PromptCStr = cell(length(CombineOverlapNumbers),1);
defAnsCStr = cell(length(CombineOverlapNumbers),1);
if(length(CombineOverlapNumbers)>=2)
    for Ind = 1:length(CombineOverlapNumbers)
        PromptCStr{Ind} = [num2str(Ind),'.Overlaps: [',JoinStr(CombineOverlapNumbers{Ind},','),']; Color:'];
        defAnsCStr{Ind} = ['[',JoinStr((Ind/length(CombineOverlapNumbers).*TemplateColor),','),']']; %['[0,0,',num2str(Ind/length(CombineOverlapNumbers)),']'];
    end
else
    PromptCStr{1} = ['1.Overlaps: [',JoinStr(CombineOverlapNumbers{1},','),']; Color:'];
    defAnsCStr{1} = '[0,0,1]';
end
AnswerColorsOverlaps = inputdlg(PromptCStr,'ColorsOverlaps',1,defAnsCStr);

ColorsOverlaps = zeros(length(AnswerColorsOverlaps),3);
for Ind = 1:length(AnswerColorsOverlaps)
    ColorsOverlaps(Ind,:) = str2num(AnswerColorsOverlaps{Ind});
end

end

%% Join Numbers as Strings with Separators
function [JoinedStr] = JoinStr(NumVector,SeparatorStr)
if(length(NumVector)>=2)
    JoinedStr = [num2str(NumVector(1)),SeparatorStr,JoinStr(NumVector(2:end),SeparatorStr)];
else
    JoinedStr = num2str(NumVector);
end

end