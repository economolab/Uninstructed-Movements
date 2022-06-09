function delPSTH = getPSTHbyDel(params,met,obj)
delPSTH.left = cell(1,length(params.delay));         % (1 x number of delay lengths)
delPSTH.right = cell(1,length(params.delay));
for g = 1:length(params.delay)                      % For each delay length...
    gix = find(met.del_trialid{1}==g);              % Get the trial IDs in the first condition that have the current delay length
    tempPSTH = mean(obj.trialpsth_cond{1}(:,:,gix),3); tempPSTH = squeeze(tempPSTH);
    delPSTH.right{g} = tempPSTH;

    gix = find(met.del_trialid{2}==g);              % Get the trial IDs in the first condition that have the current delay length
    tempPSTH = mean(obj.trialpsth_cond{2}(:,:,gix),3); tempPSTH = squeeze(tempPSTH);
    delPSTH.left{g} = tempPSTH;
end
end