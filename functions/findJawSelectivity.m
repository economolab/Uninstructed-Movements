function [sel_indivSess] = findJawSelectivity(obj,met,conditions,taxis,params)

jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'prob',params);    % (1 x conditions cell array)
% Each cell: (time x trials in that condition)
jawvel.right = mean(jaw_by_cond{1},2,'omitnan');
jawvel.left = mean(jaw_by_cond{2},2,'omitnan');

sel_indivSess = jawvel.right - jawvel.left;

% findJawSelectivity
