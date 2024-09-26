function orthmode = orthogonalizeModes(mode, kinfns)
kinmodes = zeros(numel(mode.(kinfns{1})),numel(kinfns));
for i = 1:numel(kinfns)
    kinmodes(:,i) = mode.(kinfns{i});
end
orthModes = gschmidt(kinmodes);

for i = 1:numel(kinfns)
    orthmode.(kinfns{i}) = orthModes(:,i);
end
end  % orthogonalizeModes