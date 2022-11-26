function [psth,fromsess] = concatPSTH(objs)
psth = objs{1}.psth;
% Get the number of cells from each session
[~,numCells_indiv,~] = size(objs{1}.psth);
fromsess = ones(1,numCells_indiv);
for i = 2:numel(objs)
    psth = cat(2,psth,objs{i}.psth);
    [~,numCells,~] = size(objs{i}.psth);
    session = i*ones(1,numCells);
    fromsess = cat(2,fromsess,session);
end


end % concatPSTH