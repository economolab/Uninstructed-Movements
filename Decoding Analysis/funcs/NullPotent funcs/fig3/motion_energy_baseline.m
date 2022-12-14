clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error left, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM --- 
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end


%% motion energy baseline

me_presample = [];
me_sample = [];
me_delay = [];
me_gocue = [];
medist = [];
for i = 1:numel(me)
    medist = [medist ; me(i).data(:)];

    sample = obj(i).bp.ev.sample;
    delay = obj(i).bp.ev.delay;
    gocue = obj(i).bp.ev.goCue;

    for trix = 1:size(me(i).data,2)
        tm = obj(i).time + gocue(trix);
        % presample
        ix1 = find(tm>=(sample(trix)-0.5),1,'first');
        ix2 = find(tm<=(sample(trix)),1,'last');
        dat = me(i).data(ix1:ix2,trix);
        me_presample = [me_presample ; dat];
        % sample
        ix1 = find(tm>=(sample(trix)),1,'first');
        ix2 = find(tm<=(delay(trix)),1,'last');
        dat = me(i).data(ix1:ix2,trix);
        me_sample = [me_sample ; dat];
        % delay
        ix1 = find(tm>=(delay(trix)),1,'first');
        ix2 = find(tm<=(gocue(trix)),1,'last');
        dat = me(i).data(ix1:ix2,trix);
        me_delay = [me_delay ; dat];
        % gocue
        ix1 = find(tm>=(gocue(trix)),1,'first');
        ix2 = numel(tm);
        dat = me(i).data(ix1:ix2,trix);
        me_gocue = [me_gocue ; dat];
    end

end


%%
close all
clear h



f = figure; hold on
h(1) =     histogram(medist);
h(end+1) = histogram(me_gocue);
h(end+1) = histogram(me_delay);
h(end+1) = histogram(me_sample);
h(end+1) = histogram(me_presample);


% generate colors
ncols = numel(h);
hh = [0 28 60 120 235] ./ 360;
s = linspace(1,1,ncols);
v = linspace(1,1,ncols);
for i = 1:ncols
    hsv(i,:) = [hh(i) s(i) v(i)];
end
cols = hsv2rgb(hsv);

as = [0.6 0.6 0.4 0.3 0.3 ];
for i = 1:numel(h)
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = as(i);
    h(i).Normalization = 'countdensity';
    h(i).FaceColor = cols(i,:);
    h(i).EdgeAlpha = 0.1;
end

xlim([-2,100])
f.Color = 'k';
ax = gca;
ax.Color = 'k';
ax.XAxis.Color = 'w';
ax.YAxis.Color = 'w';

xlabel('motion energy')
ylabel('counts')

ax.FontSize = 12;

ll = zeros(5, 1);
for i = 1:numel(h)
    ll(i) = plot(NaN,NaN,'-','Color',cols(i,:));
end
leg = legend(ll, 'all','gocue','delay','sample','presample');
leg.TextColor = 'w';

set(f, 'InvertHardCopy', 'off');

%%

close all
clear h

f = figure;
f.Position = [281         223        1422         765];
t = tiledlayout('flow');
toplot = {'medist','me_gocue','me_delay','me_sample','me_presample'};
for i = 1:numel(toplot)
    ax = nexttile;
    h(i) = histogram(eval(toplot{i}));
    title(toplot{i},'Interpreter','none')
end

% generate colors
ncols = numel(h);
hh = [0 28 120 235 310] ./ 360;
s = 0.75;
s = linspace(s,s,ncols);
v = linspace(1,1,ncols);
for i = 1:ncols
    hsv(i,:) = [hh(i) s(i) v(i)];
end
cols = hsv2rgb(hsv);

for i = 1:numel(h)
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 1;
    h(i).Normalization = 'countdensity';
    h(i).FaceColor = cols(i,:);
%     h(i).EdgeAlpha = 0.1;
end

% xlim([-2,100])
% f.Color = 'k';
% ax = gca;
% ax.Color = 'k';
% ax.XAxis.Color = 'w';
% ax.YAxis.Color = 'w';

xlabel(t,'motion energy')
ylabel(t,'counts')

% ll = zeros(5, 1);
% for i = 1:numel(h)
%     ll(i) = plot(NaN,NaN,'-','Color',cols(i,:));
% end
% leg = legend(ll, 'all','gocue','delay','sample','presample');
% leg.TextColor = 'w';
% 
% set(f, 'InvertHardCopy', 'off');

