function NumPossibilities = CheckPossibilities(ComponentIndex,DispResult)
% Number of possibilities for mixing below a CompnentIndex.
% Just a silly little test.
%
%Usage:
%       NumPossibilities = CheckPossibilities(ComponentIndex,DispResult);
%

%% Inputs
% ComponentIndex = 9; %1:ComponentIndex areas involved --> find the number of steps from one step below (denoted as "Lim") to the ComponentIndex.
if(length(ComponentIndex)>1) %use functional programming for recursive awesomeness!
    NumPossibilities = zeros(length(ComponentIndex),1);
    for Ind = 1:length(ComponentIndex)
        NumPossibilities(Ind) = CheckPossibilities(ComponentIndex(Ind),DispResult);
    end
    return;
end

Lim = ComponentIndex-1;

%DispResult
try
    if(isempty(DispResult))
        DispResult = 0;
    end
catch CATCH_DispResult
    disp_catch(CATCH_DispResult,[mfilename,'>CATCH_DispResult'],'CATCH_DispResult');
    DispResult = 0;
end

%% sum of choices divided by 2 should give the number of possibilities
sum = 0;
for Ind = 1:(Lim-1)
    sum = sum+nchoosek(Lim,Ind);
end
NumPossibilities = sum/2;

%% display
if(DispResult)
    disp(['Number of possibilities below ComponentIndex= ',num2str(ComponentIndex),': ',num2str(NumPossibilities)]);
end

end