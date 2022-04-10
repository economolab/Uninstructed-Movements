function featix = findDLCFeatIndex(traj,view,featName)
% traj is obj.traj
% view is 1 or 2, for side or bottom cam  respectively
% featName is a 1x1 cell containing string of feature name

% returns the index of obj.traj{view}(trialix).featNames that corresponds
% to featName input argument

[~,mask] = patternMatchCellArray(traj{view}(1).featNames,featName,'all'); % use trial 1's featNames field to find index, should be same for all trials
featix = find(mask);

if isempty(featix)
    error(['Could not find index for feature: ' featName{1}]);
end

if numel(featix) > 1
    fns = traj{view}(1).featNames;
    keep = false(size(featix));
    for i = 1:numel(featix)
        if strcmpi(featName,fns{featix(i)})
            keep(i) = true;
        end
    end
    featix = featix(keep);
end


end % findDLCFeatIndex