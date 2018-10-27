function [H,OverlapStruct] = DisplayOverlapStructData(OverlapStruct)
% This function can be used to display the overlap matrices from the OverlapStruct for inspection.
%
%
%Usage:
%       [H,OverlapStruct] = DisplayOverlapStructData(OverlapStruct);
%       [H,OverlapStruct] = DisplayOverlapStructData(); %load OverlapStruct using spm_select.
%       [H,OverlapStruct] = DisplayOverlapStructData(GenerateOverlapMatrices(ClusterAnalysisDir,'gen');); %generate OverlapStruct using GenerateOverlapMatrices.m, see help of GenerateOverlapMatrices.m.
%
%
%V1.1
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.1: (03.11.2016): some improvements in plotstyle for orientation. V1.0: (02.11.2016): initial implementation based on previous tests

%% close all figures?
if(strcmp('Yes',questdlg('Close all figures?','Close figures?','Yes','No','Yes')))
    close all
end

%% Check inputs
if(~exist('OverlapStruct','var'))
    OverlapStructPath = spm_select(1,'mat','Select OverlapStruct.mat file...');
    disp(['Loading OverlapStruct from "',OverlapStructPath,'"...']);
    load(OverlapStructPath);
else
    if(isempty(OverlapStruct))
        OverlapStructPath = spm_select(1,'mat','Select OverlapStruct.mat file...');
        disp(['Loading OverlapStruct from "',OverlapStructPath,'"...']);
        load(OverlapStructPath);
    else
        if(~isstruct(OverlapStruct))
            error('Input OverlapStruct is not a struct?! Try loading or generating OverlapStruct using GenerateOverlapMatrices.m, see help of GenerateOverlapMatrices.m');
        end
    end
end

%% get Tick Labels
[Labels.XTickVector,Labels.XTickLabels,Labels.YTickVector,Labels.YTickLabels] = GetTickLabels(OverlapStruct);
    
%% Do the display
H = {}; %init empty
H{end+1} = PlotMatrix(OverlapStruct.DiceCoeff,       Labels,[0 1],'Dice Coefficients');
H{end+1} = PlotHist(  OverlapStruct.DiceCoeff,      -1,'Hist: Dice Coefficients');

H{end+1} = PlotMatrix(OverlapStruct.MinClOverlap,    Labels,[0 1],'Minimum Cluster Overlap Coefficients');
H{end+1} = PlotHist(  OverlapStruct.MinClOverlap,   -1,'Hist: Minimum Cluster Overlap Coefficients');

H{end+1} = PlotMatrix(OverlapStruct.OverlapSize,     Labels,[   ],'Overlap Size, i.e., cluster sizes in diagonal and voxels overlapping outside of it.');
H{end+1} = PlotHist(  OverlapStruct.OverlapSize,'diag','Hist: (DIAGONAL) Overlap Size, i.e., cluster sizes in diagonal and voxels overlapping outside of it.');
H{end+1} = PlotHist(  OverlapStruct.OverlapSize,    -1,'Hist: (OFFDIAGONAL) Overlap Size, i.e., cluster sizes in diagonal and voxels overlapping outside of it.');

H{end+1} = PlotMatrix(OverlapStruct.JaccardCoeff,       Labels,[0 1],'Jaccard Coefficients');
H{end+1} = PlotHist(  OverlapStruct.JaccardCoeff,      -1,'Hist: Jaccard Coefficients');

H{end+1} = PlotMatrix(OverlapStruct.SumContribution, Labels,[0 1],'Sum Normalized StatsVals Contributions in Cluster Masks');
H{end+1} = PlotHist(  OverlapStruct.SumContribution,-1,'Hist: Sum Normalized StatsVals Contributions in Cluster Masks');

H{end+1} = PlotMatrix(OverlapStruct.AbsCorrCoeff, Labels,[0 1],'AbsCorrCoeff');
H{end+1} = PlotHist(  OverlapStruct.AbsCorrCoeff,-1,'Hist: AbsCorrCoeff');

% H{end+1} = PlotMatrix(OverlapStruct.MinClOverlap.*OverlapStruct.DiceCoeff,      Labels,[0 1],'Product MinClOverlap.*DiceCoeff');
% H{end+1} = PlotHist(  OverlapStruct.MinClOverlap.*OverlapStruct.DiceCoeff,      -1,'Hist: Product MinClOverlap.*DiceCoeff');
% 
% H{end+1} = PlotMatrix(OverlapStruct.SumContribution.*OverlapStruct.MinClOverlap,Labels,[0 1],'Product SumContribution.*MinClOverlap');
% H{end+1} = PlotHist(  OverlapStruct.SumContribution.*OverlapStruct.MinClOverlap,-1,'Hist: Product SumContribution.*MinClOverlap');
% 
% H{end+1} = PlotMatrix(OverlapStruct.JaccardCoeff.*OverlapStruct.MinClOverlap,Labels,[0 1],'Product JaccardCoeff.*MinClOverlap');
% H{end+1} = PlotHist(  OverlapStruct.JaccardCoeff.*OverlapStruct.MinClOverlap,-1,'Hist: Product JaccardCoeff.*MinClOverlap');
% 
% H{end+1} = PlotMatrix(OverlapStruct.SumContribution.*OverlapStruct.MinClOverlap.*OverlapStruct.DiceCoeff,Labels,[0 1],'Product SumContribution.*MinClOverlap.*DiceCoeff');
% H{end+1} = PlotHist(  OverlapStruct.SumContribution.*OverlapStruct.MinClOverlap.*OverlapStruct.DiceCoeff,-1,'Hist: Product SumContribution.*MinClOverlap.*DiceCoeff');

H{end+1} = PlotMatrix(OverlapStruct.SumContribution.*OverlapStruct.MinClOverlap.*OverlapStruct.DiceCoeff.*OverlapStruct.JaccardCoeff.*OverlapStruct.AbsCorrCoeff,Labels,[0 1],'Product SumContribution.*MinClOverlap.*DiceCoeff.*JaccardCoeff.*AbsCorrCoeff');
H{end+1} = PlotHist(  OverlapStruct.SumContribution.*OverlapStruct.MinClOverlap.*OverlapStruct.DiceCoeff.*OverlapStruct.JaccardCoeff.*OverlapStruct.AbsCorrCoeff,-1,'Hist: Product SumContribution.*MinClOverlap.*DiceCoeff.*JaccardCoeff.*AbsCorrCoeff');

% Do this with plotmatrix
% H{end+1} = PlotCompare(OverlapStruct.DiceCoeff,   OverlapStruct.MinClOverlap,   -1,'Dice Coefficients','Minimum Cluster Overlap Coefficients');
% H{end+1} = PlotCompare(OverlapStruct.DiceCoeff,   OverlapStruct.SumContribution,-1,'Dice Coefficients','SumNorm StatsVals Contributions');
% H{end+1} = PlotCompare(OverlapStruct.MinClOverlap,OverlapStruct.SumContribution,-1,'Minimum Cluster Overlap Coefficients','SumNorm StatsVals Contributions');
% H{end+1} = PlotCompare(OverlapStruct.DiceCoeff,   OverlapStruct.OverlapSize,    -1,'Dice Coefficients','Overlap Size');
% H{end+1} = PlotCompare(OverlapStruct.OverlapSize, OverlapStruct.SumContribution,-1,'Overlap Size','SumNorm StatsVals Contributions');
% H{end+1} = PlotCompare(OverlapStruct.MinClOverlap,OverlapStruct.OverlapSize,    -1,'Minimum Cluster Overlap Coefficients','Overlap Size');
% H{end+1} = PlotCompare(OverlapStruct.DiceCoeff,   OverlapStruct.JaccardCoeff,   -1,'Dice Coefficients','Jaccard Coefficients');
% H{end+1} = PlotCompare(OverlapStruct.MinClOverlap,OverlapStruct.JaccardCoeff,   -1,'Minimum Cluster Overlap Coefficients','JaccardCoeff');
% H{end+1} = PlotCompare(OverlapStruct.JaccardCoeff,OverlapStruct.SumContribution,-1,'JaccardCoeff','SumNorm StatsVals Contributions');
% H{end+1} = PlotCompare(OverlapStruct.JaccardCoeff,OverlapStruct.OverlapSize,    -1,'JaccardCoeff','Overlap Size');

H{end+1} = CollectDataForPlotMat(OverlapStruct);

%% Highlight cluster overlaps?
ChoiceHighlight = questdlg('Highlight Cluster Overlaps using selection criteria on Overlap Matrices?','Select Data?','Yes','No','Yes');
while(strcmp('Yes',ChoiceHighlight))
    AnswerHighlights = inputdlg({'NVoxelMin= ';'MinDiceCoeff= ';'MinMinClOverlap= ';'MinJaccardCoeff= ';'MinSumContribution= ';'MinAbsCorrCoeff= '},'Limits',1,{'10';'0.25';'0.5';'0.1';'0.25';'0.3'});
    NVoxelMin          = str2double(AnswerHighlights{1});
    MinDiceCoeff       = str2double(AnswerHighlights{2}); 
    MinMinClOverlap    = str2double(AnswerHighlights{3});
    MinJaccardCoeff    = str2double(AnswerHighlights{4});
    MinSumContribution = str2double(AnswerHighlights{5});
    MinAbsCorrCoeff       = str2double(AnswerHighlights{6});
    
    OverlapMat = (OverlapStruct.DiceCoeff>=MinDiceCoeff).*(OverlapStruct.MinClOverlap>=MinMinClOverlap).*(OverlapStruct.JaccardCoeff>=MinJaccardCoeff).*(OverlapStruct.SumContribution>=MinSumContribution).*(OverlapStruct.OverlapSize>=NVoxelMin);
    CLsizes    = diag(OverlapStruct.OverlapSize);
    Overlaps   = sum(OverlapMat,2);
    H{end+1} = PlotMatrix(OverlapMat,Labels,[0 1],['OverlapMat for NVoxelMin= ',num2str(NVoxelMin),'; MinDiceCoeff= ',num2str(MinDiceCoeff),'; MinMinClOverlap= ',num2str(MinMinClOverlap),'; MinJaccardCoeff= ',num2str(MinJaccardCoeff),'; MinSumContribution= ',num2str(MinSumContribution),'; MinAbsCorrCoeff= ',num2str(MinAbsCorrCoeff)]);
    H{end+1} = figure(); clf; 
    subplot(2,2,[1 2]); imagesc(Overlaps); title({'sum(OverlapMat,2)';['NVoxelMin= ',num2str(NVoxelMin),'; MinDiceCoeff= ',num2str(MinDiceCoeff),'; MinMinClOverlap= ',num2str(MinMinClOverlap),'; MinJaccardCoeff= ',num2str(MinJaccardCoeff),'; MinSumContribution= ',num2str(MinSumContribution),'; MinAbsCorrCoeff= ',num2str(MinAbsCorrCoeff)]}); colormap('hot'); colorbar;
    subplot(2,2,3); histogram(CLsizes(Overlaps==0),0.5+(0:1:max(CLsizes(Overlaps==0))+1)); title('CLsizes for "0" Overlaps');
    subplot(2,2,4); histogram(CLsizes(Overlaps~=0),[min(unique(CLsizes(Overlaps~=0)))-0.5;unique(CLsizes(Overlaps~=0))+0.5]); title('CLsizes for NONzero Overlaps');
    ChoiceHighlight = questdlg('Highlight Cluster Overlaps using selection criteria on Overlap Matrices?','Select Data?','Yes','No','Yes');
    if(strcmp('Yes',ChoiceHighlight))
        LastNVoxelMin          = NVoxelMin;
        LastMinDiceCoeff       = MinDiceCoeff;
        LastMinMinClOverlap    = MinMinClOverlap;
        LastMinJaccardCoeff    = MinJaccardCoeff;
        LastMinSumContribution = MinSumContribution;
        LastMinAbsCorrCoeff    = MinAbsCorrCoeff;
        LastOverlapMat = OverlapMat;
    end
end
if(exist('LastOverlapMat','var')&&strcmp('Yes',questdlg('Show diff of last two OverlapMat?','Show diff?','Yes','No','Yes')))
    H{end+1} = PlotMatrix(LastOverlapMat-OverlapMat,Labels,[-1 1],{'LastOverlapMat-OverlapMat for'; ['NVoxelMin= [',num2str(LastNVoxelMin),',',num2str(NVoxelMin),']; MinDiceCoeff= [',num2str(LastMinDiceCoeff),',',num2str(MinDiceCoeff),'];']; ['MinMinClOverlap= [',num2str(LastMinMinClOverlap),',',num2str(MinMinClOverlap),']; MinJaccardCoeff= [',num2str(LastMinJaccardCoeff),',',num2str(MinJaccardCoeff),']; MinSumContribution= [',num2str(LastMinSumContribution),',',num2str(MinSumContribution),']; MinAbsCorrCoeff= [',num2str(LastMinAbsCorrCoeff),',',num2str(MinAbsCorrCoeff),']']});
    H{end+1} = figure(); clf; imagesc(sum(LastOverlapMat,2)-sum(OverlapMat,2)); title({'sum(LastOverlapMat,2)-sum(OverlapMat,2)';['NVoxelMin= [',num2str(LastNVoxelMin),',',num2str(NVoxelMin),']; MinDiceCoeff= [',num2str(LastMinDiceCoeff),',',num2str(MinDiceCoeff),'];']; ['MinMinClOverlap= [',num2str(LastMinMinClOverlap),',',num2str(MinMinClOverlap),']; MinJaccardCoeff= [',num2str(LastMinJaccardCoeff),',',num2str(MinJaccardCoeff),']; MinSumContribution= [',num2str(LastMinSumContribution),',',num2str(MinSumContribution),']; MinAbsCorrCoeff= [',num2str(LastMinAbsCorrCoeff),',',num2str(MinAbsCorrCoeff),']']}); colormap('hot'); colorbar;
end

%% Done.
disp(' ');
disp('Done.');
disp(' ');
end

%% subfunctions
%% Plot Matrix
function H = PlotMatrix(Matrix,Labels,Limits,TitleStr)

XTickVector = Labels.XTickVector;
XTickLabels = Labels.XTickLabels;

YTickVector = Labels.YTickVector;
YTickLabels = Labels.YTickLabels;

H = figure(); clf;
if(isempty(Limits))
    imagesc(Matrix); title(TitleStr); xlabel('Clusters'); ylabel('Clusters');
else
    imagesc(Matrix,Limits); title(TitleStr); xlabel('Clusters'); ylabel('Clusters');
end
set(gca,'YTick',YTickVector); set(gca,'YTickLabel',YTickLabels); set(gca,'XTick',XTickVector); set(gca,'XTickLabel',XTickLabels); set(gca,'XTickLabelRotation',45);

end

%% PlotHist
function H = PlotHist(Data,TrilInd,TitleStr)

if(ischar(TrilInd))
    if(strcmp(TrilInd,'diag'))
        Data = diag(Data);
    end
else
    Data = Data((Data.*tril(ones(size(Data)),TrilInd))~=0);
end

NBins = round(length(Data(:))/10);
if(NBins>200)
    NBins = 200;
end

try
    H = figure(); clf;
    yyaxis left
    HistObj = histogram(Data,NBins); title(TitleStr); hold on
    
    edges  = HistObj.BinEdges;
    counts = HistObj.Values;
    counts = [counts,counts(end)];
    CumSumCounts = cumsum(counts);
    MaxCumSumCounts = max(CumSumCounts);
    Ind80Perc = find(CumSumCounts>=80.*MaxCumSumCounts/100,1);
    Ind95Perc = find(CumSumCounts>=95.*MaxCumSumCounts/100,1);
    CumSumCounts = 100.*CumSumCounts./MaxCumSumCounts;
    
    yyaxis right
    plot(edges,CumSumCounts,'rx-','LineWidth',2,'MarkerSize',6); hold on
    plot(edges(Ind80Perc).*ones(size(0:0.1:CumSumCounts(Ind80Perc))),0:0.1:CumSumCounts(Ind80Perc),'g--','LineWidth',2); hold on
    plot(edges(Ind95Perc).*ones(size(0:0.1:CumSumCounts(Ind95Perc))),0:0.1:CumSumCounts(Ind95Perc),'k--','LineWidth',2); hold on
    
    text(edges(Ind80Perc),CumSumCounts(Ind80Perc)/2,'80%','HorizontalAlignment','right');
    text(edges(Ind80Perc),CumSumCounts(Ind80Perc)/2,'20%','HorizontalAlignment','left');
    text(edges(Ind95Perc),CumSumCounts(Ind95Perc)/2,'95%','HorizontalAlignment','right');
    text(edges(Ind95Perc),CumSumCounts(Ind95Perc)/2,' 5%','HorizontalAlignment','left');
catch CATCH_PlotHist
    assignin('base','CATCH_PlotHist',CATCH_PlotHist);
    if(~exist('H','var'))
        close(H);
    end
    disp('PlotHist failed!');
    disp(CATCH_PlotHist.identifier);
    disp(CATCH_PlotHist.message);
end

end

%% PlotCompare
function H = PlotCompare(Mat1,Mat2,TrilInd,TitleStr1,TitleStr2)

if(ischar(TrilInd))
    if(strcmp(TrilInd,'diag'))
        Mat1 = diag(Mat1);
        Mat2 = diag(Mat2);
    end
else
    Selector = (((Mat1.*tril(ones(size(Mat1)),TrilInd))~=0)+((Mat2.*tril(ones(size(Mat2)),TrilInd))~=0))~=0;
    Mat1 = Mat1(Selector);
    Mat2 = Mat2(Selector);
end

H = figure(); clf;
subplot(1,2,1); 
plot(Mat1,Mat2,'kx'); hold on
plot([min([0;min(Mat1)]);max([1;max(Mat1)])],[min([0;min(Mat2)]);max([1;max(Mat2)])],'r--'); 
xlabel(TitleStr1); ylabel(TitleStr2);
title([TitleStr2,' VS ',TitleStr1]);

subplot(1,2,2); 
plot(Mat2,Mat1,'kx'); hold on
plot([min([0;min(Mat2)]);max([1;max(Mat2)])],[min([0;min(Mat1)]);max([1;max(Mat1)])],'r--');
xlabel(TitleStr2); ylabel(TitleStr1);
title([TitleStr1,' VS ',TitleStr2]);
end

%% plotmatrix(X)
function H = CollectDataForPlotMat(OverlapStruct)

Titles = {'DiceCoeff';'MinClOverlap';'JaccardCoeff';'SumContribution';'AbsCorrCoeff';'OverlapSize'};
Data = [GetData(OverlapStruct.DiceCoeff,-1),GetData(OverlapStruct.MinClOverlap,-1),GetData(OverlapStruct.JaccardCoeff,-1),GetData(OverlapStruct.SumContribution,-1),GetData(OverlapStruct.AbsCorrCoeff,-1),GetData(OverlapStruct.OverlapSize,-1)];
Selector = sum(Data~=0,2)~=0;
H = figure(); clf;
[S,AX] = plotmatrix(Data(Selector,:)); title('plotmatrix DiceCoeff MinClOverlap JaccardCoeff SumContribution AbsCorrCoeff OverlapSize');
for Ind1 = 1:size(AX,1)
    for Ind2 = 1:size(AX,2)
        if(Ind1==length(Titles))
            xlabel(AX(Ind1,Ind2),Titles{Ind2});
        end
        if(Ind2==1)
            ylabel(AX(Ind1,Ind2),Titles{Ind1});
        end
    end
end
        

end

%% GetData
function Data = GetData(Data,TrilInd)

if(ischar(TrilInd))
    if(strcmp(TrilInd,'diag'))
        Data = diag(Data);
    end
else
    Data = Data(tril(ones(size(Data)),TrilInd)~=0);
    Data = Data(:);
end

end

%% GetTickLabels
function [XTickVector,XTickLabels,YTickVector,YTickLabels] = GetTickLabels(OverlapStruct)
% Produce TickLabels and Vector

NICs     = max(OverlapStruct.Setup.ICnumClNum(:,1));
OrderICs = ceil(log10(NICs));
MaxCL    = max(OverlapStruct.Setup.ICnumClNum(:,2));
OrderCLs = ceil(log10(MaxCL));

%X
XTickVector = 1:2:size(OverlapStruct.Setup.ICnumClNum,1);
XTickLabels = cell(length(XTickVector),1);
LastICset   = 0;
for Ind = 1:length(XTickVector)
    if((OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),2)>=round(OverlapStruct.Setup.ClsPerIC(OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),1))/2))&&LastICset<OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),1))
        TickStr = ['IC ',num2str(OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),1),['%0',num2str(OrderICs),'d']),'    ',num2str(OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),2),['%0',num2str(OrderCLs),'d'])];
        LastICset = OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),1); %IC has been set in TickStr prep for the next.
    else
        TickStr = ['       ',num2str(OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),2),['%0',num2str(OrderCLs),'d'])];
    end
    XTickLabels{Ind} = TickStr;
end

%Y
YTickVector = 2:2:size(OverlapStruct.Setup.ICnumClNum,1);
YTickLabels = cell(length(YTickVector),1);
LastICset   = 0; %init 
for Ind = 1:length(YTickVector)
    if((OverlapStruct.Setup.ICnumClNum(YTickVector(Ind),2)>=round(OverlapStruct.Setup.ClsPerIC(OverlapStruct.Setup.ICnumClNum(YTickVector(Ind),1))/2))&&LastICset<OverlapStruct.Setup.ICnumClNum(YTickVector(Ind),1))
        TickStr = ['IC ',num2str(OverlapStruct.Setup.ICnumClNum(YTickVector(Ind),1),['%0',num2str(OrderICs),'d']),'    ',num2str(OverlapStruct.Setup.ICnumClNum(YTickVector(Ind),2),['%0',num2str(OrderCLs),'d'])];
        LastICset = OverlapStruct.Setup.ICnumClNum(XTickVector(Ind),1); %IC has been set in TickStr prep for the next.
    else
        TickStr = ['       ',num2str(OverlapStruct.Setup.ICnumClNum(YTickVector(Ind),2),['%0',num2str(OrderCLs),'d'])];
    end
    YTickLabels{Ind} = TickStr;
end


end
