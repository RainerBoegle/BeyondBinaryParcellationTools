function [MixModelOutDir,MMthresh,MELODICdir,MMcommStr] = CreateMixModelThreshMaps(MELODICdir,MMthresh,MixModelOutDir)
% This function calls melodic (on the command line) to do the mixture model
% fitting and thresholding, including the output of thresholded z-stats maps.
% This output can then be used for creating the Base2Encode map for display
% of all maps created from level-2 of the multi-level ICA, i.e. split of initial map.
%
%Usage: 
%      [MixModelOutDir,MMthresh,MELODICdir] = CreateMixModelThreshMaps(MELODICdir,MMthresh,MixModelOutDir);
%      [MixModelOutDir,MMthresh,MELODICdir] = CreateMixModelThreshMaps(); %select and input all necessary directories and the threshold, default is p>0.5.
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (23.08.2015): initial implementation

UseMask = 1; %default is that we use the mask if available

%% check inputs --> ask user to select dirs & threshold if not input
%CATCH_MELODICdir
try
    if(~exist(MELODICdir,'dir'))
        error(['MELODICdir: "',MELODICdir,'" not found!']);
    else
        if(~exist([MELODICdir,filesep,'melodic_IC.nii'],'file'))
            error(['"melodic_IC.nii" not found in MELODICdir: "',MELODICdir,'"!']);
        elseif(~exist([MELODICdir,filesep,'melodic_mix'],'file'))
            error(['"melodic_mix" not found in MELODICdir: "',MELODICdir,'"!']);
        elseif(~exist([MELODICdir,filesep,'mask.nii'],'file'))
            disp(['"mask.nii" not found in MELODICdir: "',MELODICdir,'"!']);
            disp( 'No problem, will do mixture model thresholding without mask, ie use implicit mask.');
            UseMask = 0;
        end
    end
catch CATCH_MELODICdir
    MELODICdir = spm_select(1,'dir','Select MELODICdir, ie the directory containing the results of melodic for fMRI data...');
    disp_catch(CATCH_MELODICdir,[mfilename,'>CATCH_MELODICdir'],'CATCH_MELODICdir'); %assignin('base','CATCH_MELODICdir',CATCH_MELODICdir);
end
%CATCH_MixModelOutDir
try
    MixModelOutDir = CheckMixModelOutDir(MixModelOutDir,MELODICdir);
catch CATCH_MixModelOutDir
    MixModelOutDir = [MELODICdir,filesep,'MixModelThreshResults'];
    MixModelOutDir = CheckMixModelOutDir(MixModelOutDir,MELODICdir);
    disp_catch(CATCH_MixModelOutDir,[mfilename,'>CATCH_MixModelOutDir'],'CATCH_MixModelOutDir'); %assignin('base','CATCH_MixModelOutDir',CATCH_MixModelOutDir);
end
%CATCH_MMthresh  
try
    if(isempty(MMthresh))
        disp('MMthresh is empty! Setting it to 0.5 (default).');
        MMthresh = 0.5;
    else
        if(~((MMthresh>0)&&(MMthresh<1))) %MMthresh has to be between 0&1!!!
            disp('MMthresh has to be between 0&1!!! Setting it to 0.5 (default).');
            MMthresh = 0.5;
        end
    end
catch CATCH_MMthresh
    disp('Setting MMthresh to 0.5 (default).');
    disp_catch(CATCH_MMthresh,[mfilename,'>CATCH_MMthresh'],'CATCH_MMthresh'); %assignin('base','CATCH_MMthresh',CATCH_MMthresh);
    MMthresh = 0.5;
end            

%% create command string for melodic mixture model (gamma-gaußian-gamma)
if(UseMask)
    %             melodic -i melodic_IC.nii -m mask.nii --ICs=melodic_IC.nii --mix=melodic_mix --Oall --report --mmthresh=0.5                   -o MixModelThreshResults
    MMcommStr = ['melodic -i melodic_IC.nii -m mask.nii --ICs=melodic_IC.nii --mix=melodic_mix --Oall --report --mmthresh=',num2str(MMthresh),' -o ',MixModelOutDir];
else
    %             melodic -i melodic_IC.nii --ICs=melodic_IC.nii --mix=melodic_mix --Oall --report --mmthresh=0.5                   -o MixModelThreshResults
    MMcommStr = ['melodic -i melodic_IC.nii --ICs=melodic_IC.nii --mix=melodic_mix --Oall --report --mmthresh=',num2str(MMthresh),' -o ',MixModelOutDir];
end


%% execute command string on command line in MELODICdir
CurrDir = pwd;
cd(MELODICdir); pause(0.1); %BugFix: needed for update
status = system(MMcommStr);
if(status~=0)
    disp('An error occured! check outputs.');
end
cd(CurrDir); pause(0.1); %BugFix: needed for update

%% Done
disp('Mixture Model estimation done.');

end

%% subfunction
function MixModelOutDir = CheckMixModelOutDir(MixModelOutDir,MELODICdir)
% This function checks if MixModelOutDir is already there or needs to be created.

if(~exist(MixModelOutDir,'dir'))
    disp(['Will create MixModelOutDir: "',MixModelOutDir,'".']);
else
    ChoiceOverwrite = questdlg(['Overwrite results in MixModelOutDir "',MixModelOutDir,'"?'],'Overwrite?','Yes','No, -new dir with current date&time.','No, -select dir.','No');
    switch(ChoiceOverwrite)
        case 'No, -select dir.'
            MixModelOutDir = spm_select(1,'dir','Select output directory for mixture model stats results...');
        case 'No, -new dir with current date&time.'
            %use MELODICdir and add a new default dir using current date & time
            MixModelOutDir = [MELODICdir,filesep,'MixModelThreshResults_',datestr(now,'yyyymmmdd_HHMM')];
        otherwise
            %'Yes, -delete old results.'
            disp(['Will overwrite old results in "',MixModelOutDir,'".']);
    end
end
end