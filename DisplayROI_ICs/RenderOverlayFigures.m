function [FigNumsForRendering,NamesForRendering,OutputDir] = RenderOverlayFigures(FigNumsForRendering,NamesForRendering,OutputDir,OnlyBitMaps)
% This function can be used for rendering figures (specifically overlays made with slover) to pdf and eps.
%
%Usage:
%      [FigNumsForRendering,NamesForRendering,OutputDir] = RenderOverlayFigures(FigNumsForRendering,NamesForRendering,OutputDir);
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (28.08.2015): initial implementation

if(~exist('OnlyBitMaps','var'))
    OnlyBitMaps = 0;
end

%% settings
dpi = 600; %150; %draft %300; %good looking but large %600; %for papers (giant)

disp('Trying to render overlays...');
%% check inputs
%FigNumsForRendering: exlude empty assignment as this can happen in a cell & allow error for this input because we need it, can't and shouldn't automatically determine the figure numbers!
tmp = []; %init empty --> needed later in case the next if block is not run
if(iscell(FigNumsForRendering))
    tmp = FigNumsForRendering; clear FigNumsForRendering
    NonEmptyInds = []; %init empty and keep track of non-empty entries in cell
    FigNumsForRendering = [];
    Ind = 1; %init
    while(Ind<=length(tmp))
        if(~isempty(tmp{Ind}))
            FigNumsForRendering = [FigNumsForRendering; tmp{Ind}];
            NonEmptyInds = [NonEmptyInds; Ind];
        end
        Ind = Ind+1; %next please!
    end
end
%NamesForRendering: allow to be empty of not input --> save as Fig#.pdf
try
    if(isempty(NamesForRendering))
        NamesForRendering = cell(length(FigNumsForRendering),1);
        disp('Creating names for figure rendering outputs from figure numbers...');
        for Ind = 1:length(FigNumsForRendering)
            NamesForRendering{Ind} = ['Fig',num2str(FigNumsForRendering(Ind)),'.pdf'];
        end
    else
        if(iscell(NamesForRendering))
            if(length(NamesForRendering)==length(tmp))
                if(~isempty(NonEmptyInds))
                    disp('Using the supplied names for figure rendering outputs that are corresponding to non-empty figure number assignments...');
                    NamesForRendering = NamesForRendering(NonEmptyInds); %assign those fitting using non-empty indices
                end
            else
                if(length(NamesForRendering)==length(FigNumsForRendering))
                    disp('Using input names for figure rendering outputs...');
                end
            end
        else
            if(ischar(NamesForRendering)) %try to convert to cellstr
                disp('Input names given as char/string. Converting to cellstr, hoping for the best...');
                NamesForRendering = cellstr(NamesForRendering);
            end
        end
    end
catch CATCH_NamesForRendering
    disp_catch(CATCH_NamesForRendering,[mfilename,'>CATCH_NamesForRendering'],'CATCH_NamesForRendering');
    NamesForRendering = cell(length(FigNumsForRendering),1);
    disp('Trying to create names for figure rendering outputs from figure numbers...');
    for Ind = 1:length(FigNumsForRendering)
        NamesForRendering{Ind} = ['Fig',num2str(FigNumsForRendering(Ind)),'.pdf'];
    end
end
%OutputDir --> allow empty or missing, -then use working directory
try
    if(isempty(OutputDir))
        disp(['Will save in current working directory ("',pwd,'")...']);
        OutputDir = pwd;
    else
        if(ischar(OutputDir))
            if(~exist(OutputDir,'dir'))
                disp(['Will create OutputDir "',OutputDir,'"...']);
                mkdir(OutputDir);
            else
                disp(['Will save outputs to "',OutputDir,'"...']);
            end
        else
            if(iscellstr(OutputDir))
                tmp2 = OutputDir{1}; clear OutputDir
                OutputDir = tmp2;
                if(~exist(OutputDir,'dir'))
                    disp(['Will create OutputDir "',OutputDir,'"(cell{1}-->char)...']);
                    mkdir(OutputDir);
                else
                    disp(['Will save outputs to "',OutputDir,'"(cell{1}-->char)...']);
                end
            else
                disp('WARNING: OutputDir is neither char or cellstr! Using current working directory instead...');
                OutputDir = pwd;
            end
        end
    end
catch CATCH_OutputDir
    disp_catch(CATCH_OutputDir,[mfilename,'>CATCH_OutputDir'],'CATCH_OutputDir');
    disp(['WARNING: Will save in current working directory ("',pwd,'")...']);
    OutputDir = pwd;
end                                

%% try to render the Overlays using opengl & painters for pdf & eps output
PCstr = computer;
if(strcmpi(PCstr(1:5),'pcwin'))
    if(OnlyBitMaps~=0)
        if(OnlyBitMaps==1)
            %tiff
            disp('Rendering tiffs...');
            try
                for IndPlot= 1:length(FigNumsForRendering)
                    print(FigNumsForRendering(IndPlot),[OutputDir,filesep,NamesForRendering{IndPlot},'.tiff'],'-dtiff','-noui','-painters',['-r',num2str(dpi)]);
                end
            catch CATCH_tiffs
                disp_catch(CATCH_tiffs,[mfilename,'>CATCH_tiffs'],'CATCH_tiffs');
            end
        end
        
        %'-dpng'
        disp('Rendering pngs...');
        try
            for IndPlot= 1:length(FigNumsForRendering)
                print(FigNumsForRendering(IndPlot),[OutputDir,filesep,NamesForRendering{IndPlot},'.png'],'-dpng','-noui','-painters',['-r',num2str(dpi)]);
            end
        catch CATCH_tiffs
            disp_catch(CATCH_tiffs,[mfilename,'>CATCH_tiffs'],'CATCH_tiffs');
        end
    else
        %pdf
        disp('Rendering pdfs...');
        try
            for IndPlot= 1:length(FigNumsForRendering)
                print(FigNumsForRendering(IndPlot),[OutputDir,filesep,NamesForRendering{IndPlot},'.pdf'],'-dpdf','-noui','-opengl',['-r',num2str(dpi)],'-loose');
            end
        catch CATCH_pdfs
            disp_catch(CATCH_pdfs,[mfilename,'>CATCH_pdfs'],'CATCH_pdfs');
        end
        
        %tiff
        disp('Rendering tiffs...');
        try
            for IndPlot= 1:length(FigNumsForRendering)
                print(FigNumsForRendering(IndPlot),[OutputDir,filesep,NamesForRendering{IndPlot},'.tiff'],'-dtiff','-noui','-painters',['-r',num2str(dpi)]);
            end
        catch CATCH_tiffs
            disp_catch(CATCH_tiffs,[mfilename,'>CATCH_tiffs'],'CATCH_tiffs');
        end
        
        %-depsc2
        disp('Rendering eps''...');
        try
            for IndPlot= 1:length(FigNumsForRendering)
                print(FigNumsForRendering(IndPlot),[OutputDir,filesep,NamesForRendering{IndPlot},'.ps'],'-depsc2 ','-noui','-painters',['-r',num2str(dpi)]);
            end
        catch CATCH_eps
            disp_catch(CATCH_eps,[mfilename,'>CATCH_eps'],'CATCH_eps');
        end
    end
else
    %ps
    disp('Rendering ps...');
    try
        for IndPlot= 1:length(FigNumsForRendering)
            print(FigNumsForRendering(IndPlot),[OutputDir,filesep,NamesForRendering{IndPlot},'.ps'],'-dpsc','-noui','-opengl',['-r',num2str(dpi)],'-loose');
        end
    catch CATCH_ps
        disp_catch(CATCH_ps,[mfilename,'>CATCH_ps'],'CATCH_ps');
    end
end

%% DONE.
disp('done.');

end