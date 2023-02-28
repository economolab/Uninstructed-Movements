% Find CDContext aligned to context switch trials
function toAW_CDCont = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, CDCont,switchtype,start,stop)
BoundaryTrix = (2*nBufferTrials)+1;              % Current switch trial will be middle, then want to look at nBufferTrials before and nBufferTrials after switch
nSwitchTrix = length(obj(sessix).(switchtype));       % Num switch trials in this session
toAW_CDCont = NaN(BoundaryTrix,nSwitchTrix);     % (# trials you want to look at before and after switches x # of switch trials)
centerix = nBufferTrials + 1;
for ii = 1:nSwitchTrix
    trix = obj(sessix).(switchtype)(ii);             % Get the current switch trial
    temp = mean(CDCont(start:stop,trix),1,'omitnan');     % Take the average presample CDContext on this switch trial
    toAW_CDCont(centerix,ii) = temp;            % Store this in toAW_CDCont variable
    for b = 1:nBufferTrials
        ix = trix+b;                            % Find the trial that is 'b' number of trials after the switch trial
        if ix > obj(sessix).bp.Ntrials          % If the trial num is greater than the number of trials in the session...
            toAW_CDCont(centerix+b,ii) = NaN;   % Take the avg presample CDContext to be NaN
        else
            temp = mean(CDCont(start:stop,ix),1,'omitnan');   % Take the average presample CDContext on this trial
            toAW_CDCont(centerix+b,ii) = temp;
        end

        clear ix temp
        ix = trix-b;                            % Find the trial that is 'b' number of trials before the switch trial
        if ix < 1                               % If the trial num is negative (aka not possible)...
            toAW_CDCont(centerix+b,ii) = NaN;   % Take the avg presample CDContext to be NaN
        else
            temp = mean(CDCont(start:stop,ix),1,'omitnan');     % Take the average presample CDContext on this trial
            toAW_CDCont(centerix-b,ii) = temp;
        end
    end
end
toAW_CDCont = toAW_CDCont';
end