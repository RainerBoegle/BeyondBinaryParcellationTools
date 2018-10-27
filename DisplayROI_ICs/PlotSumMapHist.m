function [Hfig,ExpBinCenters,Counts,PercentVox,SumMapPath] = PlotSumMapHist(SumMapPath)
% This function showes the histogram of the SumMap 
% (NB: input data are all nonzero voxels, expected bin centers are from 1 to max(unique(data))).
%
% How many voxels are singular, -only in one component will be "1", 
% All others ("2...N") indicate how many components are mixed, 
% i.e. how many voxels are a mix of two components and so on.
%
%Usage:
%       [Hfig,ExpBinCenters,Counts,PercentVox,SumMapPath] = PlotSumMapHist(SumMapPath); %SumMapPath can be missing then spm_select is used for manual selection and this selection is output in "SumMapPath".
%
%
%V1.05
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.05: (06.10.2016): small bug fix to prevent extra output of figure object contents on command window. V1.0: (28.08.2015): initial implementation

%% load SumMap NIFTI data
try
    V_SumMap = spm_vol(SumMapPath);
catch CATCH_loadSumMap
    disp_catch(CATCH_loadSumMap,[mfilename,'>CATCH_loadSumMap'],'CATCH_loadSumMap'); %couldn't load it use spm_select for manual loading
    SumMapPath = spm_select(1,'image','Select "SumMap.nii" NIFTI-file...');
    V_SumMap   = spm_vol(SumMapPath);
end

%% get filename and directory for title 
[SumMapDir,SumMapfName,ext] = fileparts(SumMapPath);
SumMapfName = [SumMapfName,ext];

%% get data
NIIdata = V_SumMap.private.dat(:);
NIIdata(NIIdata==0) = []; %remove zeros
Nnonzero = length(NIIdata);

%% set expected bin centers
ExpBinCenters = 1:max(unique(NIIdata(:)));

%% plot histogram and get counts for output using hist function (for better handling than bar function)
%% also plot the percentage covered by these kinds of voxels
Counts = hist(NIIdata(:),ExpBinCenters); %need to do it again or hist will suppress the plotting.

PercentVox = zeros(size(Counts)); %percent covered by these kinds of voxels
prev = 0; %init previous count at zero
for Ind = 1:length(Counts)
    PercentVox(Ind) = Counts(Ind)+prev;
    prev = PercentVox(Ind);
end
PercentVox = 100.*PercentVox./Nnonzero;

Htmp = figure; 
[AX,H1,H2] = plotyy(ExpBinCenters,Counts,ExpBinCenters,PercentVox,'bar','plot');
set(get(AX(1),'Ylabel'),'String','NVoxels for "m"-components mixed'); 
set(get(AX(2),'Ylabel'),'String','% covered [cumulative]'); 
set(AX(2),'YLim',[0 102]);
xlabel('Number of components mixed together ("m")');
set(H2,'Marker','x','LineStyle','-','Color','r','LineWidth',2,'MarkerSize',12);
title({['Histogram of "',SumMapfName,'" (',num2str(Nnonzero),' non-zero Voxels).']; ['From directory: "',SumMapDir,'"']},'Interpreter','none');
%add text description of individual percentages
Tmp = get(AX(1)); %need to get the axis highlighted such that text is added correctly.
for Ind = 1:length(Counts)
    if(round(100.*Counts(Ind)./Nnonzero)>=10)
        yPos  = 0.95*Counts(Ind);
    else
        if(round(100.*Counts(Ind)./Nnonzero)>.5)
            yPos = 1.1*Counts(Ind);
        else
            yPos = 0.5*Nnonzero/100; %fix at 0.5%
        end
    end
    if(round(100.*Counts(Ind)./Nnonzero)>=.1)
        Str = [num2str(round(1000.*Counts(Ind)./Nnonzero)./10),'%'];
    else
        Str = [num2str(round(100000.*Counts(Ind)./Nnonzero)./1000),'%'];
    end
    text(ExpBinCenters(Ind),yPos,Str,'HorizontalAlignment','center','EdgeColor',[1 0 0],'BackgroundColor',[1 1 1],'Color',[0 0 0]);
end

Hfig    = cell(3,1);
Hfig{1} = Htmp;
Hfig{2} = H1;
Hfig{3} = H2;


%% done.
disp('SumMapHist plot done.');

end

