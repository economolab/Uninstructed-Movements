%% Remove unwanted sessions
[meta,objs] = UseInclusionCritera(objs,meta);

%% Selectivity in jaw movements 
% Find the selectivity in prob of jaw movement for each session, separated
% by delay length
conditions = {1,2};
sel_indivSess = cell(1,length(meta));
for s = 1:length(meta)
    obj = objs{s};
    met = meta(s);
    [sel_indivSess{s}] = findHazardedJawSelectivity(obj,met,conditions,taxis,params);
end
%% Average selectivity in jaw movements first across sessions from an animal, then across animals
deltouse = 3;

[animNames,uc,nAnimals] = getAnimalNames(meta);
sessbyAnm = cell(1,nAnimals);
nConditions = numel(conditions);

sel_allAnm = NaN(length(taxis),nAnimals);
for a = 1:nAnimals
    currAnm = uc(a);    % Get the current animal name
    temp = strcmp(animNames,currAnm);     % Find the entries in 'meta' and 'obj' that correspond to this animal
    sessID = find(temp);
    nSess = length(sessID); sessbyAnm{a} = nSess;           % Number of sessions for this animal
    
    touse = NaN(length(taxis),nSess);
    for n = 1:nSess
        currsess = sessID(n);
        touse(:,n) = sel_indivSess{currsess}{deltouse};
    end
    sel_allAnm(:,a) = medfilt1(mean(touse,2,'omitnan'),10);
end

finalsel = mean(sel_allAnm,2,'omitnan');

figure();
plot(taxis,finalsel,'black','LineWidth',2)
ylabel('Selectivity in jaw prob (%)')
xlabel('Time before delay onset')
title('Jaw selectivity across all hazard del animals')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATIC DELAY SESSIONS

% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
ctrlmeta = [];
ctrlmeta = loadJEB6_ALMVideo(ctrlmeta);
ctrlmeta = loadJEB7_ALMVideo(ctrlmeta);
ctrlmeta = loadEKH1_ALMVideo(ctrlmeta);
ctrlmeta = loadEKH3_ALMVideo(ctrlmeta);
ctrlmeta = loadJGR2_ALMVideo(ctrlmeta);
ctrlmeta = loadJGR3_ALMVideo(ctrlmeta);

taxis = ctrlmeta(end).tmin:ctrlmeta(end).dt:ctrlmeta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
ctrlobjs = loadObjs(ctrlmeta);

for i = 1:numel(ctrlmeta)
    obj = ctrlobjs{i};
    obj.condition = params.condition;
    % get trials and clusters to use
    ctrlmeta(i).trialid = findTrials(obj, obj.condition);   % Get which trials pertain to the behavioral conditions you are looking at
    cluQuality = {obj.clu{ctrlmeta(i).probe}(:).quality}';  % Get clusters that are of the qualities that you specified
    ctrlmeta(i).cluid = findClusters(cluQuality, ctrlmeta(i).quality);
    % align data
    obj = alignSpikes(obj,ctrlmeta(i),params);              % Align the spike times to the event that you specified
    % get trial avg psth, single trial data, and single trial data grouped
    % by condition (aka R 2AFC, R AW, etc.)
    obj = getPSTHs(obj,ctrlmeta(i));
    ctrlobjs{i} = obj;
end

%% Remove unwanted sessions
[ctrlmeta,ctrlobjs] = UseInclusionCritera(ctrlobjs,ctrlmeta);

%% Selectivity in jaw movements 
% Find the selectivity in prob of jaw movement for each session, separated
% by delay length
ctrlsel_indivSess = cell(1,length(ctrlmeta));
for s = 1:length(ctrlmeta)
    obj = ctrlobjs{s};
    met = ctrlmeta(s);
    [ctrlsel_indivSess{s}] = findJawSelectivity(obj,met,conditions,taxis,params);
end
%% Average selectivity in jaw movements first across sessions from an animal, then across animals

[animNames,uc,nAnimals] = getAnimalNames(ctrlmeta);
sessbyAnm = cell(1,nAnimals);
nConditions = numel(conditions);

ctrlsel_allAnm = NaN(length(taxis),nAnimals);
for a = 1:nAnimals
    currAnm = uc(a);    % Get the current animal name
    temp = strcmp(animNames,currAnm);     % Find the entries in 'meta' and 'obj' that correspond to this animal
    sessID = find(temp);
    nSess = length(sessID); sessbyAnm{a} = nSess;           % Number of sessions for this animal
    
    touse = NaN(length(taxis),nSess);
    for n = 1:nSess
        currsess = sessID(n);
        touse(:,n) = ctrlsel_indivSess{currsess};
    end
    ctrlsel_allAnm(:,a) = medfilt1(mean(touse,2,'omitnan'),10);
end

ctrlfinalsel = mean(ctrlsel_allAnm,2,'omitnan');
figure();
plot(taxis,ctrlfinalsel,'black','LineWidth',2)
ylabel('Selectivity in jaw prob (%)')
xlabel('Time before delay onset')
title('Jaw selectivity across all static del animals')
%% sqrt[(R - L selectivity)^2]
finalsel = sqrt(finalsel.^2);
ctrlfinalsel = sqrt(ctrlfinalsel.^2);

figure();
plot(taxis,mySmooth(finalsel,70),'black','LineWidth',3)
hold on;
plot(taxis,mySmooth(ctrlfinalsel,70),'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',3)
xline(0,'black','LineStyle','--','LineWidth',2)
xlim([-1.3 0.9])
legend('Random Del (1.2 s)','Static Del','Location','best','FontSize',13.5)
ylabel('Squared selectivity in jaw prob (%)','FontSize',13.5)
xlabel('Time before delay onset','FontSize',13.5)
title('Avg jaw selectivity across animals','FontSize',14)