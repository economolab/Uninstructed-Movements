function delPSTH = getPSTHbyDel(params,del,obj)
delPSTH.left = cell(1,length(params.delay));         % (1 x number of delay lengths)
delPSTH.right = cell(1,length(params.delay));
for g = 1:length(params.delay)                      % For each delay length...
    gix = find(del.del_trialid{1}==g);              % Get the trial IDs in the first condition that have the current delay length
    condPSTH = obj.trialdat(:,:,params.trialid{1});
    tempPSTH = mean(condPSTH(:,:,gix),3); tempPSTH = squeeze(tempPSTH);
    delPSTH.right{g} = tempPSTH;

    gix = find(del.del_trialid{2}==g);              % Get the trial IDs in the first condition that have the current delay length
    condPSTH = obj.trialdat(:,:,params.trialid{2});
    tempPSTH = mean(condPSTH(:,:,gix),3); tempPSTH = squeeze(tempPSTH);
    delPSTH.left{g} = tempPSTH;
end
end