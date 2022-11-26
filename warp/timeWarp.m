function [obj,lick] = timeWarp(obj,params,prbnum)

% get first params.Nlicks lick times, for each trial
[lickStart,lickEnd,lickDur] = findLickTimes(obj,params.nLicks);
lick.lickStart = lickStart;
lick.lickEnd = lickEnd;
lick.lickDur = lickDur;

%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
%find median lick times for each lick across trials
% only using trials where a lick was present in median calculation
med = findMedianLickTimes(lickStart,lickEnd,lickDur, params.nLicks);
lick.med = med;
%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
% find fit for each trial and each lick
p = trialWarpFits(lickStart,lickEnd,lickDur,med,obj,params);

%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
% warp spike times for each cluster (only spike times between go cue and first params.nLicks licks get
% warped)
for cluix = 1:numel(obj.clu{params.probe(prbnum)})
    
%     disp(['Warping spikes for cluster ' num2str(cluix) ' / ' num2str(numel(obj.clu{params.probe(prbnum)}))])
    obj.clu{params.probe(prbnum)}(cluix).trialtm_warped = obj.clu{params.probe(prbnum)}(cluix).trialtm;

    for trix = 1:obj.bp.Ntrials
        % find spike times for current trial
        spkmask = ismember(obj.clu{params.probe(prbnum)}(cluix).trial,trix);
        spkix = find(spkmask);
        spktm = obj.clu{params.probe(prbnum)}(cluix).trialtm(spkmask);
        
        gc = obj.bp.ev.goCue(trix);
        ls = lickStart{trix}; % current trial lick(lix) start times
        le = lickEnd{trix}; % current trial lick(lix) end times
        ld = lickDur{trix}; % current trial lick(lix) durations
        
        
        for lix = 1:numel(lickStart{trix}) % lick index for current trial
            p_cur = p{trix,lix}; % current fit parameters for current trial and lick number
            
             % find spike ix b/w to warp
            if lix==1
                mask = (spktm>=gc) & (spktm<=le(lix));
            else
                mask = (spktm>le(lix-1)) & (spktm<=le(lix));
            end
            
            tm = spktm(mask);
            
            % warp
            warptm = polyval(p_cur,tm);
            obj.clu{params.probe(prbnum)}(cluix).trialtm_warped(spkix(mask)) = warptm;
        end
        
    end
    
end


end % timeWarp




























