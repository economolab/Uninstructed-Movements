% OUTPUTS:
% mov.earlyMoveTrial = (nTrials x 1) logical for whether each trial is an
% early move trial or not
% mov.moveTime = (nTrials x 1) cell array where each cell contains an array
% of all the times in the trial when the animal is moving 
% mov.stationaryTime = (nTrials x 1) cell array where each cell contains an array
% of all the times in the trial when the animal is stationary

% me.data = (nTrials x 1) cell array where each cell contains the
% motionEnergy values for each time point during a trial
% me.moveThresh = value above which the animal is considered to be moving

function [meta,mov,me] = assignEarlyTrials(obj,meta,params)

% load motion energy data, assumes it's stored in same location as data obj
meta.mefn = ['motionEnergy_' meta.anm '_' meta.date];
temp = load(fullfile(meta.datapth,meta.mefn));
me = temp.me;
clear temp
if ~isstruct(me)
    warning('Assign a threshold for the motion energy to label moving vs. stationary')
end

if contains(meta.anm,'JGR2')
    obj.bp.Ntrials = obj.bp.Ntrials - 1;
end

mov.earlyMoveTrial = false(obj.bp.Ntrials,1);
mov.moveTime = cell(obj.bp.Ntrials,1);
mov.stationaryTime = cell(obj.bp.Ntrials,1);
for trix = 1:obj.bp.Ntrials
    vidtime = (1:numel(me.data.data{trix}))./400;
    
    movemask = me.data.data{trix} > me.moveThresh;
    mov.moveTime{trix} = vidtime(movemask)';
    
    stationarymask = ~movemask;
    mov.stationaryTime{trix} = vidtime(stationarymask)';
    
    delayStartEnd = [obj.bp.ev.delay(trix) ; obj.bp.ev.goCue(trix)-params.dt];
    % if moving for move than 15% of delay period, label as a
    % earlyMoveTrial
    total = sum( (vidtime >= delayStartEnd(1)) & (vidtime <= delayStartEnd(2)) );
    moving = sum( (mov.moveTime{trix} >= delayStartEnd(1)) & (mov.moveTime{trix} <= delayStartEnd(2)) );
    mov.earlyMoveTrial(trix) = (moving/total) >= params.moveThresh;
end

disp(['Number of early move trials:  ' num2str(sum(mov.earlyMoveTrial)) '/' num2str(obj.bp.Ntrials) ...
    ' (' num2str(sum(mov.earlyMoveTrial)/obj.bp.Ntrials) '%)'])


end