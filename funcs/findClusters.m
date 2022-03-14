function idx = findClusters(qualityList, qualities)
% find idx where qualityList contains at least one of the patterns in
% qualities

if any(strcmpi(qualities,'all')) || strcmpi(qualities{1},'all')
    idx = 1:numel(qualityList);
    return
end

% handle unlabeled cluster qualities
for i = 1:numel(qualityList)
    if isempty(qualityList{i})
        qualityList(i) = {'nan'};
    end
    qualityList(i) = lower(qualityList(i));
end

[~,mask] = patternMatchCellArray(qualityList, qualities, 'any');

idx = find(mask);
end % findClusters