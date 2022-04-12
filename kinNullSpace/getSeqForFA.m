function obj = getSeqForFA(obj,params)

edges = params.tmin:params.dt:params.tmax;
obj.time = edges + params.dt/2;
obj.time = obj.time(1:end-1);

% get psths by condition
obj.psth = zeros(numel(obj.time),numel(params.cluid),numel(params.condition));
for i = 1:numel(params.cluid)
    curClu = params.cluid(i);
    for j = 1:numel(params.condition)
        trix = params.trialid{j};
        spkix = ismember(obj.clu{params.probe}(curClu).trial, trix);
        
        % if no spikes found for current set of trials (trix), move on to
        % next set
        if all(~spkix)
            continue
        end
        
        if params.timeWarp
            N = histc(obj.clu{params.probe}(curClu).trialtm_aligned_warped(spkix), edges);
        else
            N = histc(obj.clu{params.probe}(curClu).trialtm_aligned(spkix), edges);
        end
        N = N(1:end-1);
        
        obj.psth(:,i,j) = mySmooth(N./numel(trix)./params.dt, params.smooth);  % trial-averaged separated by trial type
        
    end
end

% get single trial data
obj.trialdat = zeros(numel(obj.time),numel(params.cluid),obj.bp.Ntrials);
for i = 1:numel(params.cluid)
    curClu = params.cluid(i);
    for j = 1:obj.bp.Ntrials
                
        spkix = ismember(obj.clu{params.probe}(curClu).trial, j);
        
        % if no spikes found for current trial (j), move on to
        % next trial
        if all(~spkix)
            continue
        end

        if params.timeWarp
            N = histc(obj.clu{params.probe}(curClu).trialtm_aligned_warped(spkix), edges);
        else
            N = histc(obj.clu{params.probe}(curClu).trialtm_aligned(spkix), edges);
        end
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialdat(:,i,j) = N./params.dt;

    end
end


end % getSeq