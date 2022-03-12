function obj = alignSpikes(obj,meta,params)

if strcmp(params.alignEvent,'firstLick')
    % get first lick time for left and right licks
    temp = obj.bp.ev.lickL;
    idx = ~cellfun('isempty',temp);
    outL = zeros(size(temp));
    outL(idx) = cellfun(@(v)v(1),temp(idx));
    temp = obj.bp.ev.lickR;
    idx = ~cellfun('isempty',temp);
    outR = zeros(size(temp));
    outR(idx) = cellfun(@(v)v(1),temp(idx));
    firstLick = zeros(size(temp));
    % firstLick = min(outL,outR), except when outL||outR == 0
    outL(outL==0) = nan;
    outR(outR==0) = nan;
    firstLick = nanmin(outL,outR);
    firstLick(isnan(firstLick)) = 0;
    obj.bp.ev.(params.alignEvent) = firstLick;
end


% align spikes to params.alignEvent
for clu = 1:numel(obj.clu{meta.probe})                       % For each cluster...   
    
    event = NaN(length(obj.clu{meta.probe}(clu).trial),1);   % Pre-allocate a vector that will allow you to put the timing of the alignEvent for the trial that each spike is found in  
    for t = 1:numel(event)                                   % For all spikes from the given cluster...
        p = obj.clu{meta.probe}(clu).trial(t);               % Set p to the trial number that the spike is found in
        if p < obj.bp.Ntrials                                % ACCOUNTS FOR DISCREPENCIES BETWEEN #TRIALS IN SLGX AND BPOD
            event(t) = obj.bp.ev.(params.alignEvent)(p);     % Find the time of the event you are aligning to (i.e. GoCue) for the given trial
        end
    end
    
    %Time within the trial that each spike occurs relative to the
    %alignEvent
    obj.clu{meta.probe}(clu).trialtm_aligned = obj.clu{meta.probe}(clu).trialtm - event;    % Subtract the timing of the alignEvent from the spike time for all spikes
end

end % alignSpikes

