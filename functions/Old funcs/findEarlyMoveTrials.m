% Identifies trials that are considered to have uninstructed early
% movements of the jaw or tongue.  Combines them into a field called
% 'earlyMoveix' within the 'obj' variable

function obj = findEarlyMoveTrials(obj)
obj = findEarlyTongue(obj);     % If tongue is visible in more than 15% of delay epoch frames
obj = findEarlyJaw(obj);        % If jaw is moving in more than 15% of delay epoch frames

% Combine early tongue and early jaw indices into 'earlyMoveix'
earlyMoveix = unique(cat(1,obj.earlyTongueix, obj.earlyJawix));
obj.earlyMoveix = earlyMoveix;
obj.pctEarlyMove = length(obj.earlyMoveix)/size(obj.trialpsth,3);
end