function [presampleFR, presampleSigma] = baselineFRForFA(obj,params)

% baselineFR - (nCells,nTrialTypes) mean FR for each clu and tt

edges = -0.5:params.dt:mode(obj.bp.ev.sample);
time = edges + params.dt/2;
time = time(1:end-1);

% get psths by condition
psth = zeros(numel(time),numel(params.cluid),numel(params.condition));
for i = 1:numel(params.cluid)
    curClu = params.cluid(i);
    for j = 1:numel(params.condition)
        trix = params.trialid{j};
        spkix = ismember(obj.clu{params.probe}(curClu).trial, trix);

        N = histc(obj.clu{params.probe}(curClu).trialtm(spkix), edges);
        N = N(1:end-1);

        psth(:,i,j) = mySmooth(N./numel(trix)./params.dt, params.smooth);  % trial-averaged separated by trial type
    end
end

presampleFR = squeeze(median(psth,1));
presampleSigma = squeeze(std(psth,[],1));

end
