function testsingleproj = getSingleTrialProjs_TTSplit(mode,touse)

%%%% Project single trials onto the coding dimension that you
%%%% specified
fns = fieldnames(touse);
for c = 1:length(fns)
    cond = fns{c};
    nTrials = size(touse.(cond),3);
    condproj = [];
    for trix = 1:nTrials                                % For each trial...
        temp = touse.(cond)(:,:,trix);                  % Get the PSTH for all cells on that trial
        condproj = [condproj,mySmooth((temp*mode),21)]; % Project the trial PSTH onto the mode that you specified
    end
    TrixProj{c} = condproj;
end
testsingleproj = TrixProj;                               % Assign this back to the rez structure
end
