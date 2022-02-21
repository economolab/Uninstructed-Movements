function jawOnset = alignJawOnset(view,feat,obj,opts)

% figure;
jawOnset = obj.bp.ev.goCue; % start with go cue, replace with jaw onset time or keep as go cue if can't find it
for trix = 1:obj.bp.Ntrials
    tm = obj.traj{view}(trix).frameTimes;
    fs = 1./median(diff(tm));
    vidx = obj.traj{view}(trix).ts(:, 1, feat);
    vidy = obj.traj{view}(trix).ts(:, 2, feat);
    
    % get jaw pos relative to origin, fill nans, mean center
    jawPos = sqrt((vidx-nanmean(vidx)).^2 + (vidy-nanmean(vidy)).^2);
    jawPos = fillmissing(jawPos, 'nearest');
    jawPos = jawPos-mean(jawPos);
    
    % filter jaw position
    Wn = opts.f_cut/fs/2;
    [b, a] = butter(opts.f_n, Wn);
    filtJawPos = filtfilt(b, a, jawPos);
    
    % find peaks in jaw position with min peak dist and min pk prominance
    [~, locs] = findpeaks(filtJawPos, 'MinPeakDistance', ceil(opts.minpkdist*fs), 'MinPeakProminence',opts.minpkprom);
    
    try
        jtms = tm(locs);
        
        goCue = obj.bp.ev.goCue(trix)+0.5;
        
        firstPeakIX = find(jtms>goCue+0.05, 1, 'first');
        
        peaktm = jtms(find(jtms>goCue+0.05, 1, 'first'));
        troughtm = tm(find(tm(1:end-1)<peaktm&diff(filtJawPos)<0, 1, 'last'));
        peakix = find(tm<peaktm, 1, 'last');
        troughix = find(tm>troughtm, 1, 'first');
        
        thresh = mean([filtJawPos(peakix),filtJawPos(troughix)]);
        crosstm = tm(find(filtJawPos<thresh&tm<peaktm, 1, 'last'));
        
        
        
%         hold on;
        
%         plot(tm-crosstm, trix*20+filtJawPos, 'r'); grid on;
        
        %             plot(tm(1:end-1)-jtms(firstPeakIX), i*30+diff(filtJawPos), 'b');
        %             plot(tm(locs), i*30+filtJawPos(locs), 'r.', 'MarkerSize', 10);
        %             plot(goCue, i*30, 'k.', 'MarkerSize', 15);
        %             plot(jtms(firstPeakIX), i*30+zeros(size(firstPeakIX)), 'b.', 'MarkerSize', 15);
        
        
        jawOnset(trix) = crosstm;
        
    catch
        continue
    end
    
end
jawOnset = jawOnset - 0.5;

end