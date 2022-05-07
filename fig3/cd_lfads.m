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
params.lfads_run = 'run3'; % 'run3' , leave empty to use most recent run
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




% % plot psths by condition as a sanity check
% figure;
% for i = 1:size(psth,2)
%     clf; hold on
%     plot(obj.time,squeeze(psth(:,i,1)),'b','LineWidth',3)
%     plot(obj.time,squeeze(psth(:,i,2)),'r','LineWidth',3)
%     pause; hold off
% end

% % plot factors_avg by condition as a sanity check
% figure;
% for i = 1:size(factors_avg,2)
%     clf; hold on
%     plot(obj.time,squeeze(factors_avg(:,i,1)),'b','LineWidth',3)
%     plot(obj.time,squeeze(factors_avg(:,i,2)),'r','LineWidth',3)
%     pause; hold off
% end

%% coding dimensions with psth


rez.time = obj.time;
rez.psth = psth;
rez.factors = factors_avg;
rez.condition = condition;
rez.alignEvent = params.alignEvent;
rez.ev = obj.bp.ev;


%% cd late mode

e1 = -0.5;
e2 = -0.01;

[~,e1] = min(abs(obj.time - e1));
[~,e2] = min(abs(obj.time - e2));

mu = squeeze(mean(psth(e1:e2,:,:),1));
sd = squeeze(std(psth(e1:e2,:,:),[],1));

cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
rez.cdLate_mode = cd./sum(abs(cd)); % (ncells,1)


%% cd go mode

e1 = 0.01;
e2 = 0.4;

[~,e1] = min(abs(obj.time - e1));
[~,e2] = min(abs(obj.time - e2));

mu = squeeze(mean(psth(e1:e2,:,:),1));
sd = squeeze(std(psth(e1:e2,:,:),[],1));

cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
rez.cdGo_mode = cd./sum(abs(cd)); % (ncells,1)

%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(params.cluid),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end


%% PLOTS

clrs = getColors();


% plot correct trials alone
plt.title = 'Correct Trials';
plt.legend = {'Right Hit','Left Hit'};
plt.conditions = [1,2];
plt.lw = [4 4];
plt.smooth = 1;
plt.colors = {clrs.rhit,clrs.lhit};
plt.save = 0;
plotAllCDs(rez, rez.ev, params.alignEvent, plt) 

% % plot correct trials and error trials
% plt.title = 'Correct and Error Trials';
% plt.legend = {'Right Hit','Left Hit','Right Error', 'Left Error'};
% plt.conditions = [1,2,5,6];
% plt.lw = [4 4 3 3];
% plt.smooth = 51;
% plt.colors = {clrs.rhit,clrs.lhit,clrs.rmiss,clrs.lmiss};
% plotAllCDs(rez, obj.bp.ev, params.alignEvent, plt) 
% 



%% selectivity


cond = [1,2];
selpsth = psth;
latent_cdlate = zeros(size(selpsth,1),numel(cond));
latent_cdgo = zeros(size(selpsth,1),numel(cond));
for i = 1:numel(cond)
    c = cond(i);
    latent_cdlate(:,i) = selpsth(:,:,c)*rez.cdLate_mode;
    latent_cdgo(:,i) = selpsth(:,:,c)*rez.cdGo_mode;
end

sm = 1;

psth_selectivity = mySmooth(selpsth(:,:,1) - selpsth(:,:,2),sm);
cdlate_selectivity = mySmooth(latent_cdlate(:,1) - latent_cdlate(:,2),sm);
cdgo_selectivity = mySmooth(latent_cdgo(:,1) - latent_cdgo(:,2),sm);


sample = mode(rez.ev.sample) - mode(rez.ev.(params.alignEvent));
delay  = mode(rez.ev.delay) - mode(rez.ev.(params.alignEvent));

lw = 4;
figure; 
plot(rez.time,mean(psth_selectivity.^2,2),'k','LineWidth',lw)
hold on
plot(rez.time,mean(cdlate_selectivity.^2,2),'g','LineWidth',lw)
plot(rez.time,mean(cdgo_selectivity.^2,2),'m','LineWidth',lw)
plot(rez.time,(cdgo_selectivity.^2 + cdlate_selectivity.^2),'c','LineWidth',lw)

xline(sample,'k--','LineWidth',0.5);
xline(delay,'k--','LineWidth',0.5);
xline(0,'k--','LineWidth',0.5);

xlabel('Time (s) from go cue')
ylabel('Squared Selectivity')
legend('Total selectivity','late','go','late + go')
xlim([rez.time(1)+0.2,rez.time(end)])
ax = gca;
ax.FontSize = 20;


%% correlation pop selectivity vector

corr_matrix_selectivity = zeros(size(psth_selectivity,1),size(psth_selectivity,1));

for i = 1:size(corr_matrix_selectivity,1)
    for j = 1:size(corr_matrix_selectivity,1)
        temp = corrcoef(psth_selectivity(i,:),psth_selectivity(j,:));
        corr_matrix_selectivity(i,j) = temp(1,2);
    end
end

%% plot selectivity correlation matrix

figure; hold on;
imagesc(obj.time,obj.time,corr_matrix_selectivity);
colorbar; caxis([0 max(max(corr_matrix_selectivity))]);

lw = 4;
xline(sample,'w--','LineWidth',lw); yline(sample,'w--','LineWidth',lw)
xline(delay,'w--','LineWidth',lw); yline(delay,'w--','LineWidth',lw)
xline(0,'w--','LineWidth',lw); yline(0,'w--','LineWidth',lw)

xlim([rez.time(1)+0.2,rez.time(end)]);
ylim([rez.time(1)+0.2,rez.time(end)])

ax = gca;
ax.FontSize = 20;
hold off
colormap(hot)



%% coding dimensions with factors

%% cd late mode

e1 = -0.5;
e2 = -0.01;

[~,e1] = min(abs(obj.time - e1));
[~,e2] = min(abs(obj.time - e2));

mu = squeeze(mean(factors_avg(e1:e2,:,:),1));
sd = squeeze(std(factors_avg(e1:e2,:,:),[],1));

cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
rez.cdLate_mode = cd./sum(abs(cd)); % (ncells,1)


%% cd go mode

e1 = 0.01;
e2 = 0.4;

[~,e1] = min(abs(obj.time - e1));
[~,e2] = min(abs(obj.time - e2));

mu = squeeze(mean(factors_avg(e1:e2,:,:),1));
sd = squeeze(std(factors_avg(e1:e2,:,:),[],1));

cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
rez.cdGo_mode = cd./sum(abs(cd)); % (ncells,1)


%% orthogonalize

clear modes
[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(rez.cdLate_mode),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end


%% PLOTS

clrs = getColors();


% plot correct trials alone
plt.title = 'Correct Trials';
plt.legend = {'Right Hit','Left Hit'};
plt.conditions = [1,2];
plt.lw = [4 4];
plt.smooth = 1;
plt.colors = {clrs.rhit,clrs.lhit};
plt.save = 0;
plotAllCDs(rez, rez.ev, params.alignEvent, plt, 'factors') 

% % plot correct trials and error trials
% plt.title = 'Correct and Error Trials';
% plt.legend = {'Right Hit','Left Hit','Right Error', 'Left Error'};
% plt.conditions = [1,2,5,6];
% plt.lw = [4 4 3 3];
% plt.smooth = 51;
% plt.colors = {clrs.rhit,clrs.lhit,clrs.rmiss,clrs.lmiss};
% plotAllCDs(rez, obj.bp.ev, params.alignEvent, plt) 
% 


%% selectivity


cond = [1,2];
selpsth = factors_avg;
latent_cdlate = zeros(size(selpsth,1),numel(cond));
latent_cdgo = zeros(size(selpsth,1),numel(cond));
for i = 1:numel(cond)
    c = cond(i);
    latent_cdlate(:,i) = selpsth(:,:,c)*rez.cdLate_mode;
    latent_cdgo(:,i) = selpsth(:,:,c)*rez.cdGo_mode;
end

sm = 1;

psth_selectivity = mySmooth(selpsth(:,:,1) - selpsth(:,:,2),sm);
cdlate_selectivity = mySmooth(latent_cdlate(:,1) - latent_cdlate(:,2),sm);
cdgo_selectivity = mySmooth(latent_cdgo(:,1) - latent_cdgo(:,2),sm);


sample = mode(rez.ev.sample) - mode(rez.ev.(params.alignEvent));
delay  = mode(rez.ev.delay) - mode(rez.ev.(params.alignEvent));

lw = 4;
figure; 
plot(rez.time,mean(psth_selectivity.^2,2),'k','LineWidth',lw)
hold on
plot(rez.time,mean(cdlate_selectivity.^2,2),'g','LineWidth',lw)
plot(rez.time,mean(cdgo_selectivity.^2,2),'m','LineWidth',lw)
plot(rez.time,(cdgo_selectivity.^2 + cdlate_selectivity.^2),'c','LineWidth',lw)

xline(sample,'k--','LineWidth',0.5);
xline(delay,'k--','LineWidth',0.5);
xline(0,'k--','LineWidth',0.5);

xlabel('Time (s) from go cue')
ylabel('Squared Selectivity')
legend('Total selectivity','late','go','late + go')
xlim([rez.time(1)+0.2,rez.time(end)])
ax = gca;
ax.FontSize = 20;


%% correlation pop selectivity vector

corr_matrix_selectivity = zeros(size(psth_selectivity,1),size(psth_selectivity,1));

for i = 1:size(corr_matrix_selectivity,1)
    for j = 1:size(corr_matrix_selectivity,1)
        temp = corrcoef(psth_selectivity(i,:),psth_selectivity(j,:));
        corr_matrix_selectivity(i,j) = temp(1,2);
    end
end

%% plot selectivity correlation matrix

figure; hold on;
imagesc(obj.time,obj.time,corr_matrix_selectivity);
colorbar; caxis([0 max(max(corr_matrix_selectivity))]);

lw = 4;
xline(sample,'w--','LineWidth',lw); yline(sample,'w--','LineWidth',lw)
xline(delay,'w--','LineWidth',lw); yline(delay,'w--','LineWidth',lw)
xline(0,'w--','LineWidth',lw); yline(0,'w--','LineWidth',lw)

xlim([rez.time(1)+0.2,rez.time(end)]);
ylim([rez.time(1)+0.2,rez.time(end)])

ax = gca;
ax.FontSize = 20;
hold off
colormap(hot)

































