function [obj,meta] = getPsthByCond(meta,obj,params)



% get psths by condition
obj.psth = zeros(numel(obj.time),numel(meta.cluNum),numel(obj.condition));
for i = 1:numel(meta.cluNum)
    curClu = meta.cluNum(i);
    for j = 1:numel(obj.condition)
        trix = meta.trialNum{j};
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);

        obj.psth(:,i,j) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);  % trial-averaged separated by trial type
    end
end

% get single trial data
nCondTrials = sum(cellfun(@numel, meta.trialNum));
obj.trialpsth = zeros(numel(obj.time),numel(meta.cluNum),nCondTrials);
obj.trialcounts = obj.trialpsth;
for i = 1:numel(meta.cluNum)
    trialct = 1;
    for j = 1:obj.bp.Ntrials
        curClu = meta.cluNum(i);
        trix = j;
        if ~(ismember(trix,meta.trialNum{1}) || ismember(trix,meta.trialNum{2}))
            continue
        end
        
        obj.trialId(trialct) = trix;
        
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialcounts(:,i,trialct) = N;
        obj.trialpsth(:,i,trialct) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);
        
        trialct = trialct + 1;

    end
end

end % getPsthByCond











