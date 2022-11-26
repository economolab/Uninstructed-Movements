function firstLickTimes = alignTraj(thisobj,alignEvent)
if strcmp(alignEvent,'firstLick')
    % get first lick time for left and right licks
    obj = thisobj;
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
    firstLickTimes = firstLick;
end
end % alignTraj