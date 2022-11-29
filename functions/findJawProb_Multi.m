function [jaw_allAnm,nAnimals,uc,sessbyAnm] = findJawProb_Multi(objs,meta,conditions,taxis,params)

[animNames,uc,nAnimals] = getAnimalNames(meta);
sessbyAnm = cell(1,nAnimals);
nConditions = numel(conditions);

jaw_allAnm.haz = cell(1,nAnimals);
for a = 1:nAnimals
    currAnm = uc(a);    % Get the current animal name
    temp = strcmp(animNames,currAnm);     % Find the entries in 'meta' and 'obj' that correspond to this animal
    sessID = find(temp);
    nSess = length(sessID); sessbyAnm{a} = nSess;           % Number of sessions for this animal
    jawvel.haz = cell(1,nSess);             % (nSessions x num delay lengths)
    for i = 1:nSess               % For each session for this animal...
        currSess = sessID(i);
        met = meta(currSess);
        obj = objs{currSess};

        % Find the probability of jaw [trident] movement at all time points in the session for trials of
        % specific conditions
        jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'prob',params);    % (1 x conditions cell array)
        % Each cell: (time x trials in that condition)
        
        jawvel.haz{i} = mean(jaw_by_cond{1},2,'omitnan');         

    end
 
    tempjaw.haz = [];

    for i = 1:nSess
        tempjaw.haz = [tempjaw.haz,jawvel.haz{i}];
    end
    meanCurrAnm.haz = mean(tempjaw.haz,2,'omitnan');
    jaw_allAnm.haz{a} = meanCurrAnm.haz;
end


end % findHazardedJaw_Multi
