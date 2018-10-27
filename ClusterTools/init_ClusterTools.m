function [] = init_ClusterTools()
% initialize clustering tools 
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (27.04.2015): initial implementation

%% get function dir & add
InitFunPath = fileparts(mfilename('fullpath'));
addpath(InitFunPath);
disp('Clustering Tools have been added to path.');

%% is disp_catch available?
DispCatchFun = cellstr(which('disp_catch','-all'));
if(isempty(DispCatchFun))
    addpath([InitFunPath,filesep,'ServiceFunctions']);
end

end