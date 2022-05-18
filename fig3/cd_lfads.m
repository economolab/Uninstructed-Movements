clear,clc,close all


addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta.anm = 'JEB7'; % 'JEB7'  'EKH3'  'JGR2'
meta.date = '2021-04-29'; % '2021-04-29'  '2021-08-11'  '2021-11-16'


% --SPECIFY PATH TO DATA HERE
if ispc
    meta.datapth = 'C:\Code\uninstructedMovements-Munib/data';
else
    meta.datapth = '/Users/Munib/Documents/Economo-Lab/data/';
end
% the data should be stored in the following structure:
% /params.datapth/
%  --- /DataObjects/
%  --- --- /animal_name_number/ (contains data_structure*.mat)
%  --- /lfads/
%  --- --- /input/ (contains anm_name_session_date*.mat/.h5)
%  --- --- /output/ (contains model_train/valid_anm_name_session_date*.h5)


% --SPECIFY PREPROCESSING PARAMETERS YOU WANT TO CHANGE BELOW THE FUNCTION CALL
% if loading lfads data, these default params will be overwritten by params
% used to generate lfads input data
params = getDefaultParams();
% params.probe = 2; % change default params.probe from '1' to '2'


% --SPECIFY METHOD OF DENOISING SINGLE TRIALS
% we need denoised, smooth single trial neural activity. We have two
% options:
% 1) load data that's ALREADY been passed through an lfads model
% 2) perform Factor Analysis on binned single trial data, followed by
%    smoothing
params.lfads_or_fa = 'lfads'; % 'lfads' or 'fa'
params.lfads_run = 'run12'; % 'run3' , leave empty to use most recent run
params.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
params.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is
assert(strcmpi(params.full_or_reduced,'reduced'),'method to use full dimensional data doesnt exist, params.full_or_reduced should be set to `reduced`')

% --SPECIFY TIME POINTS AND LAG TO USE
params.prep = [-2.5 -0.05]; % initial and final time points (seconds) defining prep epoch, relative to alignevent
params.move = [-2.5 1.5];   % initial and final time points (seconds) defining move epoch, relative to alignevent
params.advance_movement = 0.025; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

% getNeuralActivity() returns 4 main variables
% - dat: contains lfads/fa smoothed firing rates, factors, and trial numbers
%        dat.factors and dat.rates are size (time,factors/clusters,trials)
% - meta: session meta data
% - params: parameters used for preprocessing lfads input data
% - obj: preprocessed data obj
[meta,params,obj,dat] = getNeuralActivity(meta,params);


%%

clear condition
condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
trialid = findTrials(obj,condition);

mask = ismember(dat.trials,trialid{1});
trialid{1} = find(mask);
mask = ismember(dat.trials,trialid{2});
trialid{2} = find(mask);

psth = zeros(size(dat.rates,1),size(dat.rates,2),numel(trialid));
factors_avg = zeros(size(dat.factors,1),size(dat.factors,2),numel(trialid));
for i = 1:numel(trialid)
    psth(:,:,i) = mean(dat.rates(:,:,trialid{i}),3);
    factors_avg(:,:,i) = mean(dat.factors(:,:,trialid{i}),3);
end


%% coding dimensions with psth


rez.time = obj.time;
rez.data = psth;
% rez.data = factors_avg;
rez.condition = condition;
rez.alignEvent = params.alignEvent;
rez.ev = obj.bp.ev;

%% cd early mode

e1 = mode(rez.ev.sample) - mode(rez.ev.(params.alignEvent));
e2 = mode(rez.ev.sample) - mode(rez.ev.(params.alignEvent)) + 0.4;

times = rez.time>e1 & rez.time<e2;
tempdat = rez.data(:,:,[1,2]);
mu = squeeze(mean(tempdat(times,:,:),1));
sd = squeeze(std(tempdat(times,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cdEarly_mode = cd;

%% cd late mode

times = rez.time>-0.41 & rez.time<-0.01;

tempdat = rez.data(:,:,[1,2]);
mu = squeeze(mean(tempdat(times,:,:),1));
sd = squeeze(std(tempdat(times,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cdLate_mode = cd;


%% cd go mode

times = rez.time>0.01 & rez.time<0.41;

tempdat = rez.data(:,:,[1,2]);
mu = squeeze(mean(tempdat(times,:,:),1));
sd = squeeze(std(tempdat(times,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cdGo_mode = cd;

%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(size(rez.data,2),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end

%% projections and normalize
% when pooling trajectories across sessions from hidehikos ppn paper:
% CD_late projections normalized by mean activity just before go cue
% (-0.1<t<t_go)
% CD_go projections normalized by mean activity after go cue
% (t_go<t<0.4)

normTimes{1} = rez.time>e1 & rez.time<e2; % sample
normTimes{2} = rez.time>-0.6 & rez.time<0; % delay
normTimes{3} = rez.time>0 & rez.time<0.4; % go

cond = [1 2];
for i = 1:numel(fns)
    tempmode = rez.(fns{i});
    for j = 1:numel(cond)
        c = cond(j);
        
        tempdat = rez.data(:,:,c)*rez.(fns{i});
        
        normfactor = abs(nanmean(tempdat(normTimes{i})));
        normfactor = 1;
        
        rez.([fns{i}(1:end-5) '_latent'])(:,j) = tempdat ./ normfactor;
    end
end

clear cond

%% variance explained

for i = 1:numel(fns)
    psth = rez.data;
    datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
    datacov(isnan(datacov)) = 0;
    eigsum = sum(eig(datacov));
    rez.varexp.(fns{i}(1:end-5)) = var_proj(rez.(fns{i}), datacov, eigsum);
end


%% PLOT CDs

close all
clrs = getColors();
lw = 6;
alph = 0.5;

sample = mode(rez.ev.sample - rez.ev.(params.alignEvent));
delay = mode(rez.ev.delay - rez.ev.(params.alignEvent));

sav = 0;
for i = 1:numel(fns)
    f(i) = figure; ax = axes(f(i)); hold on
    tempmean = eval(['rez.' fns{i}(1:end-5) '_latent']);
    plot(rez.time,tempmean(:,1),'Color',clrs.rhit,'LineWidth',lw);
    plot(rez.time,tempmean(:,2),'Color',clrs.lhit,'LineWidth',lw);
    
    
    xlim([rez.time(15);rez.time(end)])
    title(fns{i},'Interpreter','none')
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax = gca;
    ax.FontSize = 40;
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/cd';
        fn = [fns{i} '_JEB7_2021-04-29_run12_sm_' num2str(params.smooth)];
        mysavefig(f(i),pth,fn);
    end
    
end



%% selectivity

selpsth = rez.data;

latent_cdearly = rez.cdEarly_latent;
latent_cdlate = rez.cdLate_latent;
latent_cdgo = rez.cdGo_latent;

sm = 31;

psth_selectivity = mySmooth(selpsth(:,:,1) - selpsth(:,:,2),sm);
cdearly_selectivity = mySmooth(latent_cdearly(:,1) - latent_cdearly(:,2),sm);
cdlate_selectivity = mySmooth(latent_cdlate(:,1) - latent_cdlate(:,2),sm);
cdgo_selectivity = mySmooth(latent_cdgo(:,1) - latent_cdgo(:,2),sm);

%%
close all
sample = mode(rez(1).ev.sample) - mode(rez(1).ev.(params.alignEvent));
delay  = mode(rez(1).ev.delay) - mode(rez(1).ev.(params.alignEvent));

lw = 4;
f = figure; 
plot(rez(1).time,sum(psth_selectivity.^2,2),'k','LineWidth',lw)
hold on
plot(rez(1).time,sum(cdearly_selectivity.^2,2),'b','LineWidth',lw)
plot(rez(1).time,sum(cdlate_selectivity.^2,2),'g','LineWidth',lw)
plot(rez(1).time,sum(cdgo_selectivity.^2,2),'m','LineWidth',lw)
plot(rez(1).time,(cdearly_selectivity.^2 + cdgo_selectivity.^2 + cdlate_selectivity.^2),'c','LineWidth',lw)

xline(sample,'k--','LineWidth',0.5);
xline(delay,'k--','LineWidth',0.5);
xline(0,'k--','LineWidth',0.5);

xlabel('Time (s) from go cue')
ylabel('Squared Sum Selectivity')
legend('Total selectivity','early','late','go','early + late + go')
xlim([rez(1).time(1)+0.2,rez(1).time(end)])
ax = gca;
ax.FontSize = 20;

% create smaller axes in top right, and plot on it
axes('Position',[.2 .6 .25 .25])
box on
times = rez.time>-2.2 & rez.time<-1.2;
h2 = plot(rez(1).time(times),sum(psth_selectivity(times,:).^2,2),'k','LineWidth',lw);
hold on
plot(rez(1).time(times),sum(cdearly_selectivity(times).^2,2),'b','LineWidth',lw)
plot(rez(1).time(times),sum(cdlate_selectivity(times).^2,2),'g','LineWidth',lw)
plot(rez(1).time(times),sum(cdgo_selectivity(times).^2,2),'m','LineWidth',lw)
plot(rez(1).time(times),(cdearly_selectivity(times).^2 + cdgo_selectivity(times).^2 + cdlate_selectivity(times).^2),'c','LineWidth',lw)
% h2 = cdfplot(corrs);
% h2.LineWidth = 2;
ax = h2.Parent;
ax.XLim = [-2.2 -1.2];
% h2.Parent.XLabel.String = '';
% h2.Parent.YLabel.String = '';
% h2.Parent.Title.String = 'CDF';
h2.Parent.FontSize = 15;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/selectivity';
% fn = 'sqelectivity_JEB7_2021-04-29_run12';
% mysavefig(f,pth,fn);



%% correlation pop selectivity vector

corr_matrix_selectivity = zeros(size(psth_selectivity,1),size(psth_selectivity,1));

for i = 1:size(corr_matrix_selectivity,1)
    for j = 1:size(corr_matrix_selectivity,1)
        temp = corrcoef(psth_selectivity(i,:),psth_selectivity(j,:));
        corr_matrix_selectivity(i,j) = temp(1,2);
    end
end

%% plot selectivity correlation matrix
close all
f = figure; hold on;
imagesc(obj.time,obj.time,corr_matrix_selectivity);
colorbar; caxis([0 max(max(corr_matrix_selectivity))]);

lw = 4;
ls = '--';
col = [88, 245, 112] ./ 255;
xline(sample,ls,'Color',col,'LineWidth',lw); yline(sample,ls,'Color',col,'LineWidth',lw)
xline(delay,ls,'Color',col,'LineWidth',lw); yline(delay,ls,'Color',col,'LineWidth',lw)
xline(0,ls,'Color',col,'LineWidth',lw); yline(0,ls,'Color',col,'LineWidth',lw)

xlim([rez(1).time(1)+0.2,rez(1).time(end)]);
ylim([rez(1).time(1)+0.2,rez(1).time(end)])
xlabel('Time (s) from go cue')
ylabel('Time (s) from go cue')
ax = gca;
hold off
colormap(hot)
a = colorbar;
a.Label.String = 'Correlation of population selectivity vector';
ax.FontSize = 25;

% % 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/selectivity';
% fn = 'popselectivitycorr_pop_jeb7_2021-04-29_run12';
% mysavefig(f,pth,fn);



































