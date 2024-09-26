function obj = removeEarlyTrials_PSTH(obj,meta,params)
for c = 1:numel(obj.condition)      % For every condition...               
    trix = meta.trialid{c};         % Get the trials for that condition
    % Exclude early movement trials
    if isfield(obj,'earlyMoveix') && strcmp(params.earlytrials,'DeepLabCut')        
        cond_earlyMove = ismember(trix,obj.earlyMoveix);    % Find trials in this condition that were identified as early move
        trid = trix(~cond_earlyMove);                       % Remove those trials
        earl = trix(cond_earlyMove);
    elseif isfield(obj,'earlyMoveix') && strcmp(params.earlytrials,'motionEnergy')
        temp = find(obj.earlyMoveix);
        cond_earlyMove = ismember(trix,temp);    % Find trials in this condition that were identified as early move
        trid = trix(~cond_earlyMove);                       % Remove those trials
        earl = trix(cond_earlyMove);

    end

    obj.trialpsth_noEarly{c} = obj.trialpsth(:,:,trid); % Take only the single trial PSTHs from trials without early movements       
    obj.trialpsth_Early{c} = obj.trialpsth(:,:,earl); % Take only the single trial PSTHs from trials with early movements 
end
end