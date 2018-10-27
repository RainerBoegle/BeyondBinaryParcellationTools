function [SelectedMapsPaths,OverlapMatStruct] = CheckOutMapsBasedOnOverlapMat(OverlapMatStruct,ICoI)
% This function can be used to display IC maps based on the Overlap Matrix, i.e. based on the
% overlaps of a selected IC with others.
%
%Usage:
%       [SelectedMapsPaths,OverlapMatStruct] = CheckOutMapsBasedOnOverlapMat(OverlapMatStruct,ICoI);
%        SelectedMapsPaths                   = CheckOutMapsBasedOnOverlapMat(OverlapMatStruct,  13); %This will select the maps overlapping IC 13.       
%       [SelectedMapsPaths,OverlapMatStruct] = CheckOutMapsBasedOnOverlapMat(CreateOverlapMat(ThreshMapsList),13); %Create Overlap Matrix & then select the maps overlapping IC 13.       
%       [SelectedMapsPaths,OverlapMatStruct] = CheckOutMapsBasedOnOverlapMat(CreateOverlapMat(),13); %Create Overlap Matrix (manual selection) & then select the maps overlapping IC 13.       
% 
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (06.03.2016): initial implementation

%% get data out
OverlapMat    = OverlapMatStruct.OverlapMat; %The matrix
ICnums        = OverlapMatStruct.ICnums; %IC numbers
NVoxelsTotal  = OverlapMatStruct.NVoxelsTotal; %Total number of voxels
ThreshMapsList= OverlapMatStruct.ThreshMapsList; %List of input maps for OverlapMatrix creation

%% display overlaps for a certain ICs?
if(~exist('ICoI','var')) %select manually
    [ICoI,ok] = listdlg('ListString',cellstr(num2str(ICnums(:),'%02d')),'SelectionMode','single','ListSize',[160 min([length(ICnums)*15; 600])],'Name','IC of interest','PromptString','Select IC of interest: ','CancelString','Quit');
    if(~ok)
        SelectedMapsPaths = [];
        return;
    end
else
    if(length(ICoI)>1) %for more than one do recursive call
        SelectedMapsPaths = cell(length(ICoI),1);
        for Ind = 1:length(ICoI)
            SelectedMapsPaths{Ind} = CheckOutMapsBasedOnOverlapMat(OverlapMatStruct,ICoI(Ind));
        end
        return;
    end        
end
disp(['Displaying Overlaps with IC ',num2str(ICoI),'.']);

%% display all maps for a certain ICs?
NVoxTotStr = ['(#VoxTotal=',num2str(NVoxelsTotal),')'];

[s,sInds] = sort(OverlapMat(ICoI,:),'descend'); %sort accordint to overlaps
figure; bar(1:length(s),s); title(['(sorted) Counts for ICoI #',num2str(ICoI)]); ylabel(['#Voxels ',NVoxTotStr]);  set(gca,'XTick',1:length(s),'XTickLabel',cellstr(num2str(sInds(:))));

%% ask user including size
AllNonZeroInds = find(s~=0);
AllNonZero     = s(AllNonZeroInds);
SelStr         = cell(length(AllNonZero),1);
for Ind = 1:length(AllNonZero)
    SelStr{Ind} = [num2str(ICnums(sInds(AllNonZeroInds(Ind))),['%0',num2str(max([ceil(log10(max(ICnums(:)))); 2])),'d']),'|nVox=',num2str(OverlapMat(ICoI,sInds(AllNonZeroInds(Ind))))];
end
[SelOverlapInds,ok] = listdlg('ListString',SelStr,'ListSize',[160 min([length(ICnums(sInds(AllNonZeroInds)))*15; 600])],'Name','Overlap ICs','PromptString','Select Overlap ICs: ','CancelString','Quit');
if(~ok)
    SelectedMapsPaths = [];
    return;
end
SelOverlapICs     = ICnums(sInds(AllNonZeroInds(SelOverlapInds))); %transform this into ICnums
SelectedMapsPaths = ThreshMapsList(SelOverlapICs); %The selected paths

%% Display
DisplayStats('Stats',ThreshMapsList(ICoI)); %IC of interest
for Ind = 1:length(SelectedMapsPaths)
    DisplayStats('Stats',SelectedMapsPaths{Ind}); %IC of interest
end

%% Done.
H = helpdlg(['Done with showing overlapping ICs for IC ',num2str(ICoI,['%0',num2str(max([ceil(log10(ICoI)); 2])),'d']),'.'],'Overlap display done.');
uiwait(H);

end
