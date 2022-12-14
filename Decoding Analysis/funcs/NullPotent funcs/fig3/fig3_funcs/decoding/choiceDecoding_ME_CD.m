
%% organize data


% single trial null/potent projections
np = {'null','potent'};
cds = {'cdEarly','cdLate','cdGo'};
newcds = {'early','late','go'};
for i = 1:numel(np)
    for j = 1:numel(cds)
        cd.(np{i}).(newcds{j}) = rez.cd.(np{i}).([cds{j} '_latent']);
    end
end


% motion energy
ix = find(ismember(kin.featLeg,'motion_energy'));
motion = kinfeats{1}(:,:,ix);

% all data is (time,trials)

%% choice decoding

k = 10; % number of iterations (bootstrap)

binSize = 20; % ms
binSize = floor(binSize / (params.dt*1000)); % samples

mdlTime = obj.time(1:binSize:numel(obj.time));
numT = numel(mdlTime);




train = 0.8; % fraction of trials to use for training (1-train for testing)
cond2use = [2 3];

% n/p projections
for i = 1:numel(np)
    for j = 1:numel(cds)
        indat = cd.(np{i}).(newcds{j});
        indat(isnan(indat)) = 0;
        acc.(np{i}).(newcds{j}) = choiceDecoder(indat,numT,k,params,train,binSize,cond2use); % (time,bootiters)
    end
end
% motion eenrgy
indat = motion;
indat(isnan(indat)) = 0;
acc.me = choiceDecoder(indat,numT,k,params,train,binSize,cond2use);

%%

close all

f = figure;
f.Position = [680   110   623   868];
ax(1) = subplot(2,1,1); % null
ax(2) = subplot(2,1,2); % potent
sgtitle('Choice Decoding from CDs and motion energy')

for i = 1:numel(np)
    axes(ax(i)); hold on;
    ylabel('Accuracy')
    title(np{i})
    for j = 1:numel(newcds)
        plot(obj.time,mySmooth(mean(acc.(np{i}).(newcds{j}),2),21) ,'LineWidth',2)
    end
    plot(obj.time,mean(acc.me,2) ,'LineWidth',2,'Color','k')
end
xlabel('Time (s) from go cue')
for i = 1:numel(ax)
    ax(i).FontSize = 15;
    ax(i).XLim = [-2 2];
    ax(i).YLim = [0.4 1];
end
























































