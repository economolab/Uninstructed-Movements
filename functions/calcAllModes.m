% Function to find and orthogonalize all 7 activity modes defined in Nuo Li's preprint

function allModes = calcAllModes(obj,met,rez,params)

% stimulus mode
% 1. stimulus mode: defined during stimulus (sample) period
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'sample';
allModes.stimulus_mode = stimulusMode(obj,met,cond,epoch,rez.alignEvent);
clear cond

% choice mode
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'delay';
allModes.choice_mode = choiceMode(obj,met,cond,epoch,rez.alignEvent);
clear cond

% action mode
% 3. action mode: defined during mvmt init (0.05s before firstLick to 0.1s after)
%       (hitR - hitL) / sqrt(sum(sd for each tt ^2));

cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'action';
allModes.action_mode = actionMode(obj,met,cond,epoch,rez.alignEvent);
clear cond

% outcome mode
% 4. outcome mode: defined during response epoch (0.05s before firstLick to 1.3 s after)
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));

cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'outcome';
allModes.outcome_mode = outcomeMode(obj,met,cond,epoch,rez.alignEvent);
clear cond

% ramping mode
% 5. ramping mode: in correct trials
%       (hit_presample - hit_delay) / sqrt(sum(sd for each tt ^2));
cond{1} = params.modecondition{5};
epoch = {'presample','delay'};
allModes.ramping_mode = rampingMode(obj,met,cond,epoch,rez.alignEvent);
clear cond

% go mode
% 6. go mode: 0.1 sec before or after firstLick
%       (hit_after - hit_before) / sqrt(sum(sd for each tt ^2));
cond{1} = params.modecondition{5};
epoch = {'postgo','prego'};
allModes.go_mode = goMode(obj,met,cond,epoch,rez.alignEvent);
clear cond

% response mode
% 7. response mode:
%    a. find eigenvectors of basline subtracted PSTHs using SVD
%       aa. matrix was of size (n x (2t)), where left and right trials concatenated
%       in time
%    b. response mode = eigenvector with largest eigenvalue
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
psthcond = [1,2];
epoch = 'presample'; % used to estimate baseline firing rate
allModes.response_mode = responseMode(obj,met,cond,epoch,rez.alignEvent,psthcond);
clear cond

% orthogonalize MODES
[fns,~] = patternMatchCellArray(fieldnames(allModes),{'mode'},'all');
modes = zeros(size(obj.psth,2),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = allModes.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    allModes.(fns{i}) = orthModes(:,i);
end
end