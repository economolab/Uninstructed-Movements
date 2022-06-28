clear,clc,close all

addpath(genpath(pwd))

% finds cd early, late, go as defined in economo 2018

%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                                % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 0.02;

% smooth with causal gaussian kernel
params.smooth = 0;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'all'}; 

%% SET METADATA

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
meta = loadJEB7_ALMVideo(meta,datapth); % done

params.probe = [meta.probe]; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

use = [true false];

meta = meta(use);
params.probe = params.probe(use);

%% LOAD AND PROCESS DATA

objs = loadObjs(meta);


for metaix = 1:numel(meta)
    obj = objs{metaix};
    disp('______________________________________________________')
    disp(['Processing data for session ' [meta(metaix).anm '_' meta(metaix).date]])
    disp(' ')
    [sessparams{metaix},sessobj{metaix}] = processData(obj,params,params.probe(metaix));
end

% clean up sessparams and sessobj
for metaix = 1:numel(meta)
    params.trialid{metaix} = sessparams{metaix}.trialid;
    params.cluid{metaix} = sessparams{metaix}.cluid{params.probe(metaix)};
    
    objs{metaix} = sessobj{metaix};
    objs{metaix}.psth = objs{metaix}.psth{params.probe(metaix)};
    objs{metaix}.trialdat = objs{metaix}.trialdat{params.probe(metaix)};
    objs{metaix}.presampleFR = objs{metaix}.presampleFR{params.probe(metaix)};
    objs{metaix}.presampleSigma = objs{metaix}.presampleSigma{params.probe(metaix)};
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')

% Remove unwanted sessions

% remove sessions if:
% 1) less than 40 trials of rhit and lhit each (same as
% 2) atleast 10 clusters of quality that is not noise
use = false(size(objs));
for i = 1:numel(use)
    check(1) = numel(params.trialid{i}{1}) > 40;
    check(2) = numel(params.trialid{i}{2}) > 40;
    check(3) = numel(params.cluid{i}) >= 10;
    if all(check)
        use(i) = true;
    end
end

meta = meta(use);
objs = objs(use);
params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);



%%

params.trialid = params.trialid{1};
params.cluid = params.cluid{1};
obj = objs{1};


%% firing rates

cluix = randsample(1:numel(params.cluid),1);

nTrials = 3;
trix = randsample(1:obj.bp.Ntrials,nTrials);

%%
close all
f = figure; hold on
% f.Position = [-1379        -143        1003         914];
for i = 1:nTrials
    col = hsv2rgb([i*30,100,150]./255);
    temp = obj.clu{1}(params.cluid(cluix));
    ix = ismember(temp.trial,trix(i));
    plot(temp.trialtm(ix),i/20,'.','Color',col,'MarkerSize',10)
    
    temp2 = histc(temp.trialtm(ix),[-0.5:0.001:10]);
    temp2 = mySmooth(temp2,301);
    plot([-0.5:0.001:10],temp2+(i/19.95),'Color',col,'LineWidth',3)
end
ylim([0.025,0.18])
xlim([-0.25,5])
ax = gca;
ax.XTick = [];
ax.YTick = [];
xlabel('Time')
ylabel('Neurons')
ax.FontSize = 20;


%% one cluster, all left and right trials

cluix = randsample(1:numel(params.cluid),1);

% set conditions to calculate PSTHs for
clear cond
cond(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
cond(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off

trials = findTrials(obj,cond);

cols = getColors();
clrs{1} = cols.rhit;
clrs{2} = cols.lhit;

close all
f = figure; hold on
f.Position = [-1379        -143        1003         914];
for i = 1:numel(cond)
    subplot(2,2,i)
    temp = obj.clu{1}(params.cluid(cluix));
    ix = ismember(temp.trial,trials{i});
    plot(temp.trialtm(ix) - mode(obj.bp.ev.(params.alignEvent)),temp.trial(ix),'.','Color',clrs{i},'MarkerSize',10)
    xlim([-2.5,2.5])
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    xlabel('Time')
    ylabel('Trials')
    ax.FontSize = 20;
    
    subplot(2,2,i+2)
    
    plot(obj.time,obj.psth(:,cluix,i+1),'Color',clrs{i},'LineWidth',3)
    xlim([-2.5,2.5])
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    xlabel('Time')
    ylabel('spks / sec')
    ax.FontSize = 20;
    
    
end

%% activity in neural state space

pop = cat(1,obj.psth(:,:,2),obj.psth(:,:,3));

[~,latents] = pca(pop,'NumComponents',3);

latents = reshape(mySmooth(latents,31),numel(obj.time),2,3);
latents = permute(latents,[1 3 2]);

[~,gcix] = min(abs(obj.time - 0));

ix = 20:size(latents,1);

close all
figure; hold on
plot3(latents(ix,1,1),latents(ix,2,1),latents(ix,3,1),'Color',clrs{1}, 'LineWidth',3)
plot3(latents(ix,1,2),latents(ix,2,2),latents(ix,3,2),'Color',clrs{2}, 'LineWidth',3)

plot3(latents(gcix,1,1),latents(gcix,2,1),latents(gcix,3,1),'k.','MarkerSize',50)
plot3(latents(gcix,1,2),latents(gcix,2,2),latents(gcix,3,2),'k.','MarkerSize',50)

plot3(latents(ix(1),1,1),latents(ix(1),2,1),latents(ix(1),3,1),'.','Color',clrs{1},'MarkerSize',50)
plot3(latents(ix(1),1,2),latents(ix(1),2,2),latents(ix(1),3,2),'.','Color',clrs{2},'MarkerSize',50)

grid on;

ax = gca;
xticklabels("")
yticklabels("")
zticklabels("")
% Set major Grid lines
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = 0.5;
grid on;

xlabel('$$Dim_1$$','Interpreter','latex','FontWeight','bold')
ylabel('$$Dim_2$$','Interpreter','latex','FontWeight','bold')
zlabel('$$Dim_3$$','Interpreter','latex','FontWeight','bold')

ax.FontSize = 30;

grid on


%% pcs var explained


[~, ~, ~, ~, ve] = pca(pop);

%%
figure;
bar(ve,'EdgeColor','none')
ax = gca;
xlabel('PC')
ylabel('% VE')
ax.FontSize = 30;

% create smaller axes in top right, and plot on it
axes('Position',[.45 .4 .3 .35])
box on
h2 = cdfplot(ve);
h2.LineWidth = 4;
ax = h2.Parent;
h2.Parent.XLabel.String = '# PCs';
h2.Parent.YLabel.String = 'Cumulative %VE';
h2.Parent.Title.String = '';
h2.Parent.FontSize = 25;


figure;
plot(1:numel(ve),ve,'.-','MarkerSize',20,'LineWidth',3)
ax = gca;
xlabel('PC')
ylabel('% VE')
ax.FontSize = 30;




%% plain old pca

t = 1000;

mu1 = 5;
std1 = sqrt(500);

mu2 = 5;
std2 = sqrt(3);

n1 = mu1.*randn(t,1) + std1;
n1 = n1 - mean(n1);
n2 = (n1*1.5) + (mu2.*randn(t,1) + std2);
n2 = n2 - mean(n2);

[pcs,~,~,~,ve] = pca([n1 n2]);
pcs = gschmidt(pcs);

close all
figure; hold on;
plot(n1,n2,'.')
q1 = quiver(0,0,pcs(1,1),pcs(2,1),10);
q2 = quiver(0,0,pcs(1,2),pcs(2,2),4);
% q = quiver([0; 0],[0; 0],pcs(:,1),pcs(:,2),10);
q1.ShowArrowHead = 'off';
q1.LineWidth = 4;
q2.ShowArrowHead = 'off';
q2.LineWidth = 4;
axis equal
xlim([-20,20])
ax = gca;
xlabel('N1')
ylabel('N2')
ax.FontSize = 30;

figure;
bar(ve,'EdgeColor','none')
ax = gca;
xlabel('PC')
ylabel('% VE')
ax.FontSize = 25;

%% same but uncorrelated data


t = 1000;

mu1 = 5;
std1 = sqrt(500);

mu2 = 5;
std2 = sqrt(500);

n1 = mu1.*rand(t,1) + std1;
n1 = n1 - mean(n1);
n2 = (mu2.*rand(t,1) + std2);
n2 = n2 - mean(n2);

[pcs,~,~,~,ve] = pca([n1 n2]);
pcs = gschmidt(pcs);

close all
figure; hold on;
plot(n1,n2,'.')
q1 = quiver(0,0,pcs(1,1),pcs(2,1),2);
q2 = quiver(0,0,pcs(1,2),pcs(2,2),2);
% q = quiver([0; 0],[0; 0],pcs(:,1),pcs(:,2),10);
q1.ShowArrowHead = 'off';
q1.LineWidth = 4;
q2.ShowArrowHead = 'off';
q2.LineWidth = 4;
axis equal
% xlim([-2,2])
ax = gca;
xlabel('N1')
ylabel('N2')
ax.FontSize = 30;

figure;
bar(ve,'EdgeColor','none')
ax = gca;
xlabel('PC')
ylabel('% VE')
ax.FontSize = 25;
ylim([0 100])










