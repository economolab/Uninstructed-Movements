function jaw_allAnm = jawProb_ALLAnimals(objs,meta,taxis, conditions)

[animNames,uc,nAnimals] = getAnimalNames(meta);
nConditions = numel(conditions);

jaw_allAnm = cell(nAnimals,nConditions);

for a = 1:nAnimals      % For each animal...
    currAnm = uc(a);    % Get the current animal name
    temp = strcmp(animNames,currAnm);     % Find the entries in 'meta' and 'obj' that correspond to this animal
    sessID = find(temp);                  
    nSess = length(sessID);               % Number of sessions for this animal
    jaw_currAnm = cell(1,nSess);
    for i = 1:nSess               % For each session for this animal...
        currSess = sessID(i);   
        met = meta(currSess);
        obj = objs{currSess};
        [jaw_currAnm{i},~] = jawProbSessionAvg(obj,met,conditions);
    end
    
     
    for c = 1:nConditions
        tempjaw = [];
        for i = 1:nSess
            tempjaw = [tempjaw,jaw_currAnm{i}{c}];
        end
        meanCurrAnm = mean(tempjaw,2,'omitnan');
        jaw_allAnm{a,c} = meanCurrAnm;
    end
end
end
