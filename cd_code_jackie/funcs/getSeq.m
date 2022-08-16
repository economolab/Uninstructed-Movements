function obj = getSeq(obj,params,prbnum)

edges = params.tmin:params.dt:params.tmax;
obj.time = edges + params.dt/2;
obj.time = obj.time(1:end-1);

% get psths by condition
obj.psth{prbnum} = zeros(numel(obj.time),numel(params.cluid{prbnum}),numel(params.condition));
for i = 1:numel(params.cluid{prbnum})
    curClu = params.cluid{prbnum}(i);
    for j = 1:numel(params.condition)
        trix = params.trialid{j};
        
        spkix = ismember(obj.clu{prbnum}(curClu).trial, trix);

        % if no spikes found for current set of trials (trix), move on to
        % next set
        if all(~spkix)
            continue
        end
        
        if params.timeWarp
            N = histc(obj.clu{prbnum}(curClu).trialtm_aligned_warped(spkix), edges);
        else
            N = histc(obj.clu{prbnum}(curClu).trialtm_aligned(spkix), edges);
        end
        N = N(1:end-1);
        
        obj.psth{prbnum}(:,i,j) = mySmooth(N./numel(trix)./params.dt, params.smooth);  % trial-averaged separated by trial type
        
    end
end

% get single trial data
obj.trialdat{prbnum} = zeros(numel(obj.time),numel(params.cluid{prbnum}),obj.bp.Ntrials);
for i = 1:numel(params.cluid{prbnum})
    curClu = params.cluid{prbnum}(i);
    for j = 1:obj.bp.Ntrials
                
        spkix = ismember(obj.clu{prbnum}(curClu).trial, j);
        
        % if no spikes found for current trial (j), move on to
        % next trial
        if all(~spkix)
            continue
        end

        if params.timeWarp
            N = histc(obj.clu{prbnum}(curClu).trialtm_aligned_warped(spkix), edges);
        else
            N = histc(obj.clu{prbnum}(curClu).trialtm_aligned(spkix), edges);
        end
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialdat{prbnum}(:,i,j) = mySmooth(N./params.dt,params.smooth);

    end
end


end % getSeq