function [lickStart,lickEnd,lickDur,l_ix] = findLickTimes(obj,nLicks)

% inputs:
%   obj -> data obj
%   nLicks -> # of post go cue licks to find times for


view = 1; % use side cam
c = 1; % use x coord
feat = 1; % tongue in side cam

dt = 1/400;

% there are small discontinuities during licks (these are nans in ts)
% We can try to fill these nans by finding discontinuities within a moving
% window size of 'movWindowSize'. Then let matlab's fillmissing function do
% the work. By looking at all trials for a single session, a value of 0.025
% seconds seems to the trick. There are some licks that are not completely
% filled in and some that are merged with another, but these errors will
% likely average out over all trials
movWindowSize = 0.025 / dt; % num time points to fill discontinuities during a lick


l_ix.start = cell(obj.bp.Ntrials, 1);
l_ix.end = cell(obj.bp.Ntrials, 1);
l_ix.dur = cell(obj.bp.Ntrials, 1);

lickStart = cell(obj.bp.Ntrials, 1);
lickEnd = cell(obj.bp.Ntrials, 1); 
lickDur = cell(obj.bp.Ntrials, 1); 


%%
% figure % uncomment here and bottom of for loop to do a check
for trix = 1:obj.bp.Ntrials
    
    ts = obj.traj{view}(trix).ts(:,c,feat); % get data
    % fill discontinuities within licks (majorrity of nans between licks
    % remain nans)
    ts = fillmissing(ts,'movmedian',movWindowSize);
    t = obj.traj{view}(trix).frameTimes; % trial time
    
    % if the first time point contains a visible tongue, then we don't
    % catch this as a full lick below, so set first time point to a nan
    ts(1) = nan;
    
    gocue = obj.bp.ev.goCue(trix); % time of go cue
    
    ix_start = find(isnan(ts(1:end-1)) & ~isnan(ts(2:end)) & t(2:end)>gocue);     % ix where a lick starts (after the Go Cue)
    ix_end = zeros(size(ix_start));                                             % ix where a lick ends (after the Go cue)
    
    % assign lick end times
    % Make sure that a lick doesn't start right at the end of the trial (aka
    % wouldn't be able to find the end of the lick)
    for lix = 1:numel(ix_start)   %for each lick in the trial
        tmp = find(~isnan(ts(ix_start(lix):end-1)) & isnan(ts(ix_start(lix)+1:end)), 1, 'first')+ix_start(lix)-1;
        if isempty(tmp)
            ix_start = ix_start(1:end-1);
            ix_end = ix_end(1:end-1);
        else
            ix_end(lix) = tmp;
        end
    end
    
    % keep first nLicks licks after go cue
    if numel(ix_start) > nLicks
        ix_start = ix_start(1:nLicks);
        ix_end = ix_end(1:nLicks);
    end
    
    l_ix.start{trix} = ix_start;
    l_ix.end{trix}   = ix_end;
    l_ix.dur{trix}   = ix_end - ix_start;
    
    lickStart{trix} = ix_start .* dt;
    lickEnd{trix} = ix_end .* dt;
    lickDur{trix} = lickEnd{trix} - lickStart{trix};
    
%     plot(t,x); hold on
%     plot(t(ix1{trix}),x(ix1{trix}+1),'g.','MarkerSize',15)
%     plot(t(ix2{trix}),x(ix2{trix}-1),'r.','MarkerSize',15)
%     xlim([t(1) t(end)])
%     hold off
%     pause
    
    
end


end