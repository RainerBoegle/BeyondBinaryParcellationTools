function [] = DispOverlapMatrix(OverlapMatStruct,ICoI) %OLD: [] = DispOverlapMatrix(OverlapMat,ICnums,ICoI,NVoxelsTotal)
% This function can be used to display the overlap matrix.
%
%Usage:
%       DispOverlapMatrix(OverlapMatStruct,ICoI);
%       DispOverlapMatrix(OverlapMatStruct);
%       DispOverlapMatrix(OverlapMatStruct,13);   %plot 13th input for itself as barplot
%       DispOverlapMatrix(OverlapMatStruct,1:13); %plot 1st to 13th input for themself as barplot
%       DispOverlapMatrix(OverlapMatStruct,Inf);  %plot all inputs for themself as barplot
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (20.02.2016): initial implementation

%% NB: OverlapMatStruct
%OverlapMatStruct.
%                .OverlapMat     = OverlapMat;
%                .ThreshMapsList = ThreshMapsList;
%                .ThreshP        = ThreshP;
%                .ICnums         = ICnums;
%                .NVoxelsTotal   = NVoxelsTotal;
%                .Data.Vols        = Vols;
%                .Data.DataDim     = DataDim;
%                .Data.Mask        = Mask;
%                .Data.DataThresh2D= DataThresh2D;
%                .Data.Data2D      = Data2D;

%% get data out
OverlapMat  = OverlapMatStruct.OverlapMat; %The matrix
ICnums      = OverlapMatStruct.ICnums; %IC numbers
NVoxelsTotal= OverlapMatStruct.NVoxelsTotal; %Total number of voxels

%% figure numbers
if(size(OverlapMat,1)>41)
    bf1 = 4242;
    bf1b= 4243;
    bf2 = 8181;
else
    bf1 = 42;
    bf2 = 81;
end

%% extra info
NVoxTotStr = ['(#VoxTotal=',num2str(NVoxelsTotal),')'];

%% basic plot
figure(bf1); imagesc(OverlapMat); title(['Overlap Matrix ',NVoxTotStr]);
if(~isempty(ICnums))
    set(gca, 'XTick', 1:length(ICnums), 'XTickLabel', cellstr(num2str(ICnums(:))));
    xlabel('IC#');
else
    xlabel('Inputs');
end

figure(bf1b); imagesc(OverlapMat./(repmat(diag(OverlapMat),1,length(diag(OverlapMat))))); title(['Overlap Matrix[NORM-DIAG] ',NVoxTotStr]);
if(~isempty(ICnums))
    set(gca, 'XTick', 1:length(ICnums), 'XTickLabel', cellstr(num2str(ICnums(:))));
    xlabel('IC#');
else
    xlabel('Inputs');
end

figure(bf2); bar(diag(OverlapMat)); title('Size of each Input in voxels'); ylabel(['#Voxels ',NVoxTotStr]);

%% extra plot
if(exist('ICoI','var'))
    if(any(isinf(ICoI))) %i.e. Inf == take all
        ICoI = 1:size(OverlapMat,1);
    end
    for Ind = 1:length(ICoI)
        figure; 
        subplot(1,2,1); bar(OverlapMat(ICoI(Ind),:)); title(['Counts for ICoI #',num2str(ICoI(Ind))]); ylabel(['#Voxels ',NVoxTotStr]);
        [s,sInds] = sort(OverlapMat(ICoI(Ind),:),'descend');
        subplot(1,2,2); bar(1:length(s),s); title(['(sorted) Counts for ICoI #',num2str(ICoI(Ind))]); ylabel(['#Voxels ',NVoxTotStr]);  set(gca,'XTick',1:length(s),'XTickLabel',cellstr(num2str(sInds(:))));
    end
end

end