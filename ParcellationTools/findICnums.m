function [ok,ICnums] = findICnums(ThreshMapsList)
% Find IC nums assuming end of filenames are numbers.
% Outputs are a boolean "ok" indicating if the process was successful, i.e. "ICnums" can be trusted.
% "ICnums" are the IC numbers that have been found.
%
% These can be used for sorting the ThreshMapsList using 
% [ICnums,sortInds] = sort(ICnums);
% ThreshMapsList    = ThreshMapsList(sortInds); %InputFiles must be a cellstr otherwise this doesn't work therefore the hassle above when checking the inputs.
%
%Usage:
%      [ok,ICnums] = findICnums(ThreshMapsList); 
%
%
%V0.5
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V0.5: (22.03.2016): initial implementation WITHOUT SAFETY CHECKS FOR INPUTS


ok = 1; %assume everything is okay until we have some conflicting evidence.
ICnums = zeros(length(ThreshMapsList),1);
for Ind = 1:length(ThreshMapsList)
    [BaseDir,FName,ext] = fileparts(ThreshMapsList{Ind});
    [s,e,tok,match,tokStr,exprN,splitStr] = regexp(FName(end-4:end),'\d*'); clear s e tok tokStr exprN splitStr %keep only "match" NB: this ensures backward compatability by avoiding use of "~"
    if(isempty(match))
        ok = 0; %error!
        disp(['No IC num for Input',num2str(Ind),'! abort...']);
    else
        if(iscell(match))
            if(length(match)>1)
                ok = 0; %error!
                disp(['unclear IC num for Input',num2str(Ind),'! abort...']);
            else
                if(ischar(match{1}))
                    try
                        ICnums(Ind) = str2num(match{1});
                    catch
                        disp(['IC num conversion error found for Input',num2str(Ind),'! abort...']);
                        ok = 0; %error!
                    end
                else
                    disp(['IC num type error found for Input',num2str(Ind),'! abort...']);
                    ok = 0; %error!
                end
            end
        else
            if(ischar(match))
                try
                    ICnums(Ind) = str2num(match);
                catch
                    disp(['IC num conversion error (from char) found for Input',num2str(Ind),'! abort...']);
                    ok = 0; %error!
                end
            else
                disp(['IC num type error (from 2nd-try char) found for Input',num2str(Ind),'! abort...']);
                ok = 0; %error!
            end
        end
    end
    if(~ok)
        break;
    end
end%for-loop

end