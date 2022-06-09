function CD = orthogModes_Multi(CD, multiDelPSTH)
[fns,~] = patternMatchCellArray(fieldnames(CD),{'mode'},'all');
    modes = zeros(size(multiDelPSTH.left{1},2),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = CD.(fns{i});
    end

    orthModes = gschmidt(modes);

    for i = 1:numel(fns)
        CD.(fns{i}) = orthModes(:,i);
    end
end