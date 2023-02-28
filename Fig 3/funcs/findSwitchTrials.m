function [toAW_ix, toAFC_ix] = findSwitchTrials(bp)
% 2AFC to AW switch trials
AWix = find(bp.autowater&~bp.stim.enable);      % Find all AW trials that were not stim trials
for ii = 1:length(AWix)                         % For all AW trials...
    trix = AWix(ii);
    if trix==1 || trix==bp.Ntrials              % If the current trial is the first in the session or the last in the session
        AWix(ii) = NaN;                         % Get rid of this trial
    elseif bp.autowater(trix-1)==1              % If the previous trial was also AW
        AWix(ii) = NaN;                         % Get rid of this trial
    end
end
temp = find(isnan(AWix));
AWix(temp) = [];
toAW_ix = AWix;

% AW to 2AFC switch trials
AFCix = find(~bp.autowater&~bp.stim.enable);      % Find all 2AFC trials that were not stim trials
for ii = 1:length(AFCix)                           % For all 2AFC trials...
    trix = AFCix(ii);
    
    if trix==1 || trix==bp.Ntrials              % If the current trial is the first in the session or the last in the session
        AFCix(ii) = NaN;
    elseif bp.autowater(trix-1)==0                  % If the previous trial was also 2AFC
        AFCix(ii) = NaN;                         % Get rid of this trial
    end
end
temp = find(isnan(AFCix));
AFCix(temp) = [];
toAFC_ix = AFCix;
end
