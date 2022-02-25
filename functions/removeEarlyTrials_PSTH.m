function obj = removeEarlyTrials_PSTH(obj,meta)
for c = 1:numel(obj.condition)      % For every condition...               
    trix = meta.trialid{c};         % Get the trials for that condition
    % Exclude early movement trials
    if isfield(obj,'earlyMoveix')       
        cond_earlyMove = ismember(trix,obj.earlyMoveix);    % Find trials in this condition that were identified as early move
        trix = trix(~cond_earlyMove);                       % Remove those trials
    elseif isfield(obj,'earlyMoveTrial')
        cond_earlyMove = ismember(trix,obj.earlyMoveTrial);    % Find trials in this condition that were identified as early move
        trix = trix(~cond_earlyMove); 
    end

    obj.trialpsth_noEarly{c} = obj.trialpsth(:,:,trix); % Take only the single trial PSTHs from trials without early movements       
end
end