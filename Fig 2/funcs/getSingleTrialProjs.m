function rez = getSingleTrialProjs(rez,obj,cd)

cd2use = find(strcmp(rez(1).cd_labels,cd));             % Specify which CD you want to project onto via the input 'cd' which is a string (i.e. 'early' or 'late')

for sessix = 1:numel(obj)                               % For every session...
    %%%% Project single trials onto the coding dimension that you
    %%%% specified
    nTrials = size(obj(sessix).trialdat,3);
    TrixProj = NaN(length(rez(sessix).time),nTrials);   % time x nTrials
    mode = rez(sessix).cd_mode_orth(:,cd2use); 
    for trix = 1:nTrials                                % For each trial...
        temp = obj(sessix).trialdat(:,:,trix);          % Get the PSTH for all cells on that trial
        TrixProj(:,trix) = mySmooth((temp*mode),21);    % Project the trial PSTH onto the mode that you specified
    end

    rez(sessix).singleProj = TrixProj;                % Assign this back to the rez structure
end

end