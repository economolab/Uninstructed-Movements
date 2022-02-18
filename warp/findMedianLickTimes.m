function med = findMedianLickTimes(lickStart,lickEnd,lickDur, nLicks)

med.lickStart = nan(nLicks,1);
med.lickEnd = nan(nLicks,1);
med.lickDur = nan(nLicks,1);
for i = 1:nLicks
    % find trials with atleast (i-1) licks
    trix_w_licks = cell2mat(cellfun(@(x) (numel(x)>(i-1)), lickStart, 'UniformOutput', false)); 
    % get i_th lick start, end, and duration for each of trix_w_licks
    ith_lickStart = cell2mat(cellfun(@(x) x(i), lickStart(trix_w_licks), 'UniformOutput', false));
    ith_lickEnd = cell2mat(cellfun(@(x) x(i), lickEnd(trix_w_licks), 'UniformOutput', false));
    ith_lickDur = cell2mat(cellfun(@(x) x(i), lickDur(trix_w_licks), 'UniformOutput', false));
    % calculate median time for i_th lick
    med.lickStart(i) = median(ith_lickStart);
    med.lickEnd(i) = median(ith_lickEnd);
    med.lickDur(i) = median(ith_lickDur);
end

end
