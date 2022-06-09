function [jaw_allAnm,nAnimals,uc,sessbyAnm] = findJawProb_Multi(objs,meta,conditions,taxis,params)

[animNames,uc,nAnimals] = getAnimalNames(meta);
sessbyAnm = cell(1,nAnimals);
nConditions = numel(conditions);

jaw_allAnm.right = cell(1,nAnimals);
jaw_allAnm.left = cell(1,nAnimals);
for a = 1:nAnimals
    currAnm = uc(a);    % Get the current animal name
    temp = strcmp(animNames,currAnm);     % Find the entries in 'meta' and 'obj' that correspond to this animal
    sessID = find(temp);
    nSess = length(sessID); sessbyAnm{a} = nSess;           % Number of sessions for this animal
    jawvel.left = cell(1,nSess);         % (nSessions x number of delay lengths)
    jawvel.right = cell(1,nSess);        % (nSessions x num delay lengths)
    for i = 1:nSess               % For each session for this animal...
        currSess = sessID(i);
        met = meta(currSess);
        obj = objs{currSess};

        % Find the probability of jaw [trident] movement at all time points in the session for trials of
        % specific conditions
        jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'prob',params);    % (1 x conditions cell array)
        % Each cell: (time x trials in that condition)
        
        jawvel.right{i} = mean(jaw_by_cond{1},2,'omitnan');         
        jawvel.left{i} = mean(jaw_by_cond{2},2,'omitnan');
    end
 
    tempjaw.right = [];
    tempjaw.left = [];
    for i = 1:nSess
        tempjaw.right = [tempjaw.right,jawvel.right{i}];
        tempjaw.left = [tempjaw.left,jawvel.left{i}];
    end
    meanCurrAnm.right = mean(tempjaw.right,2,'omitnan');
    meanCurrAnm.left = mean(tempjaw.left,2,'omitnan');
    jaw_allAnm.right{a} = meanCurrAnm.right;
    jaw_allAnm.left{a} = meanCurrAnm.left;
end


end % findHazardedJaw_Multi
