function rez = getMove_SingleTrialProjs(rez,obj,cd,kin)

cd2use = find(strcmp(rez(1).cd_labels,cd));             % Specify which CD you want to project onto via the input 'cd' which is a string (i.e. 'early' or 'late')

for sessix = 1:numel(obj)                               % For every session...
    %%%% Project single trials onto the MovementCD that you specified
    nTrials = size(kin(sessix).dat,2);
    TrixProj = NaN(length(rez(sessix).time),nTrials);   % time x nTrials
    mode = rez(sessix).movecd_mode_orth(:,cd2use); 
    for trix = 1:nTrials                                % For each trial...
        temp = squeeze(kin(sessix).dat(:,trix,:));               % Get the kinematic data for all movement features on that trial
        temp = fillmissing(temp,'constant',0);
        TrixProj(:,trix) = mySmooth((temp*mode),31);    % Project the trial kinematic data onto the mode that you specified
    end

    rez(sessix).singleMoveCDProj = TrixProj;                % Assign this back to the rez structure
end

end