function CD = orthogModes(CD, obj)
[fns,~] = patternMatchCellArray(fieldnames(CD),{'mode'},'all');
    modes = zeros(size(obj.psth,2),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = CD.(fns{i});
    end

    orthModes = gschmidt(modes);

    for i = 1:numel(fns)
        CD.(fns{i}) = orthModes(:,i);
    end
end