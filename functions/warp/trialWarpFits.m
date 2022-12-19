function p = trialWarpFits(lickStart,lickEnd,lickDur,med,obj,params)


p = cell(obj.bp.Ntrials,params.nLicks);
for trix = 1:obj.bp.Ntrials
    ls = lickStart{trix}; % current trial lick start times
    le = lickEnd{trix}; % current trial lick end times
    ld = lickDur{trix}; % current trial lick durations
    
    for lix = 1:numel(ls) % lick index
        % median values for current lick
        mls = med.lickStart(lix); % median lick start time
        mle = med.lickEnd(lix); % median lick end time
        mld = med.lickDur(lix); % median lick duration
        
        % warp lick times
        % if first lick: warp from go cue to first lick end
        % if lick2 or lick>2: warp from end of previous lick to end of
        % current lick
        if lix == 1
            x = [obj.bp.ev.goCue(trix) le(lix)];
            y = [mode(obj.bp.ev.goCue) mle];
        else
            x = [le(lix-1) le(lix)];
            y = [med.lickEnd(lix-1) mle];
        end
        
        % fit original time to warped time
        p{trix,lix} = polyfit(x,y,1);
        
    end
    
end

end