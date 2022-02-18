function idx = findClusters(qualityList, qualities)
% find idx where qualityList contains at least one of the patterns in
% qualities

% handle unlabeled cluster qualities
for i = 1:numel(qualityList)
    if isempty(qualityList{i})          %If any clusters are unlabeled, put a NaN as the label
        qualityList(i) = {'nan'};
    end
end

[~,mask] = patternMatchCellArray(qualityList, qualities, 'any');    %Make a logical array that tells you which clusters are of one of the qualities that you are looking for

idx = find(mask);                       %Find the cluster numbers that are of the qualities that you are looking for
end % findClusters

