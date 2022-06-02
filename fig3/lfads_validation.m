clear,clc,close all

addpath(genpath(pwd))

% figures created in this script:
% -lfads psths vs gaussian smoothed psths (all trials and left vs. right hits)
% -visualize factors on left and right hits
% -neural reconstruction
% -correlation between lfads and gaussian smoothed psths

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta.anm = 'JEB7';          % 'JEB7'        'JGR2'         'JEB6'
meta.date = '2021-04-29';   % '2021-04-29'  '2021-11-17'   '2021-04-18'


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
params.lfads_run = 'run1'; % 'run3' , leave empty to use most recent run
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


%% lfads and gaussian smoothed single trials (all trials used to train lfads)

lfads = dat.rates;

gauss = obj.trialdat(:,:,dat.trials);

sm = 31; sigma = 30; winsize = sm*params.dt/2;
for i = 1:size(gauss,2) % for each clu
    gauss(:,i,:) = mySmooth(squeeze(gauss(:,i,:)),sm,sigma); % causal gaussian filter
end
%% save some random cells
sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

f = figure;
% ax = axes(f);
for i = 1:size(gauss,2)
    clf(f); ax = axes(f); hold on;
    
    % gaussian smoothing
    temp = squeeze(gauss(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(numel(dat.trials));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',[106, 115, 171]./255,'LineWidth',3},0.5, ax);
    
    % lfads
    temp = squeeze(lfads(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(numel(dat.trials));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',[130, 176, 120]./255,'LineWidth',3},0.5, ax);

    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
    xlabel('Time (s) from go cue')
    ylabel('Firing Rate (spikes/s)')
    title(['Cell ' num2str(params.cluid(i))])
%     legend({'LFADS','Smoothing'},'Location','best')
    ax = gca;
    ax.FontSize = 20;
    xlim([obj.time(10) obj.time(end)])
    
    pause
    hold off;
%     
%     if rand > 0.95
%         pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/lfads_psth';
%         fn = [meta.anm '_' meta.date '_cell' num2str(params.cluid(i)) 'hits_lfads_' params.lfads_run '_gausssm_' num2str(sm) '_dt_005'];
% %         fn = ['cell' num2str(params.cluid(i)) '_jeb7-2021-04-29_lfadsrun_12_sm_31_dt_0025'  ];
%         mysavefig(f,pth,fn)
%     end
    
end

%% lfads and gaussian smoothed single trials (left and right hits only)

clear temp

lfads = dat.rates;

gauss = obj.trialdat;

rhits = params.trialid{2};
lhits = params.trialid{3};

temp{1} = lfads(:,:,rhits);
temp{2} = lfads(:,:,lhits);

lfads = temp; clear temp

temp{1} = gauss(:,:,rhits);
temp{2} = gauss(:,:,lhits);

gauss = temp; clear temp

sm = 31; sigma = 10;
for j = 1:2
    for i = 1:size(gauss{j},2) % for each clu
        gauss{j}(:,i,:) = mySmooth(squeeze(gauss{j}(:,i,:)),sm,sigma); % causal gaussian filter
    end
end
%% save same cells as above
sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

clrs = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840],...
        [184, 146, 79]./255,[78, 186, 96]./255};

f = figure;
alph = 1;
for i = 1:numel(params.cluid)
    clf(f); ax = axes(f); hold on;
    
    % gaussian smoothing
    temp = squeeze(gauss{1}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(gauss{1},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{1},'LineWidth',1.5},alph, ax);
    temp = squeeze(gauss{2}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(gauss{2},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{2},'LineWidth',1.5},alph, ax);
    
    % lfads
    temp = squeeze(lfads{1}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(lfads{1},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{3},'LineWidth',3},alph, ax);
    temp = squeeze(lfads{2}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(lfads{2},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{4},'LineWidth',3},alph, ax);

    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
    xlabel('Time (s) from go cue')
    ylabel('Firing Rate (spikes/s)')
    title(['Cell ' num2str(params.cluid(i))])
%     legend({'LFADS','Smoothing'},'Location','best')
    ax = gca;
    ax.FontSize = 20;
    xlim([obj.time(10) obj.time(end)])
    
    pause
    hold off;
%     
%     if any(ismember([22,27,63],params.cluid(i)))
%         pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/lfads_psth';
%         fn = [meta.anm '_' meta.date '_cell' num2str(params.cluid(i)) 'rlhits_lfads_' params.lfads_run '_gausssm_' num2str(sm) '_dt_005' ];
% %         fn = ['cell' num2str(params.cluid(i)) '_jeb7-2021-04-29_lfadsrun_12_sm_31_dt_0025'  ];
%         mysavefig(f,pth,fn)
%     end
    
end


%% plot all cells in a big plot
clear temp

clrs = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840],...
        [184, 146, 79]./255,[78, 186, 96]./255};
    
sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));


f = figure;
alph = 1;
for i = 1:numel(params.cluid)
    ax = nexttile; hold on;
    
    % gaussian smoothing
    temp = squeeze(gauss{1}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(gauss{1},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{1},'LineWidth',1.5},alph, ax);
    temp = squeeze(gauss{2}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(gauss{2},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{2},'LineWidth',1.5},alph, ax);
    
    % lfads
    temp = squeeze(lfads{1}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(lfads{1},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{3},'LineWidth',3},alph, ax);
    temp = squeeze(lfads{2}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(lfads{2},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs{4},'LineWidth',3},alph, ax);

    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
%     xlabel('Time (s) from go cue')
%     ylabel('Firing Rate (spikes/s)')
    title(['Cell ' num2str(params.cluid(i))])
%     legend({'LFADS','Smoothing'},'Location','best')
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.FontSize = 20;
    xlim([obj.time(10) obj.time(end)])
    

%     
%     if any(ismember([22,27,63],params.cluid(i)))
%         pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/lfads_psth';
%         fn = [meta.anm '_' meta.date 'allcells'  '_rlhits_lfads_' params.lfads_run '_gausssm_' num2str(sm) '_dt_005' ];
%         mysavefig(f,pth,fn)
%     end
    
end



%% correlation b/w gauss psth and lfads psth


lfads = dat.rates;

gauss = obj.trialdat;

sm = 51; sigma = 51;
for i = 1:size(gauss,2) % for each clu
    gauss(:,i,:) = mySmooth(squeeze(gauss(:,i,:)),sm,sigma); % causal gaussian filter
end


corrs = zeros(size(gauss,2),1);
for i = 1:size(gauss,2)


    lfadspsth = mean(lfads(:,i,:),3);
    poppsth = mySmooth(mean(gauss(:,i,:),3),sm,sigma);
    
    temp = corrcoef(lfadspsth,poppsth);
    corrs(i) = temp(1,2);

end

%%

close all
f = figure; 
nbins = round(size(gauss,2) / 3);
h = histogram(corrs,nbins,'EdgeColor','none');
hold on;
xline(mean(corrs),'k','LineWidth',3)
xlabel('Correlation')
ylabel('Count')
ax = gca;
ax.FontSize = 20;

% create smaller axes in top right, and plot on it
axes('Position',[.2 .6 .3 .25])
box on
h2 = cdfplot(corrs);
h2.LineWidth = 2;
ax = h2.Parent;
h2.Parent.XLabel.String = '';
h2.Parent.YLabel.String = '';
h2.Parent.Title.String = 'CDF';
h2.Parent.FontSize = 15;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/lfads_psth_corr';
% fn = 'lfads_psthcorr_JEB7_2021-04-29_run12_gausssm_31_rhit_lhit';
% mysavefig(f,pth,fn)



%% neural reconstruction

% To compare our ability to reconstruct single-trial responses using Gaussian smoothing and LFADS, 
% we first computed peri-event time histograms (PETHs) within condition using all training trials (excluding one test trial).
% We then computed the correlation between the firing rates of each test trial (smoothed with a Gaussian kernel or 
% reconstructed with LFADS) with the PETH of the corresponding condition averaged across the training trails 
% (Figure 2—figure supplement 1A). We repeated this procedure with a different trial left out for each condition.
% We report the difference in correlation coefficient obtained after LFADS processing and Gaussian smoothing (Figure 2—figure supplement 1B).


lfads = dat.rates;

gauss = obj.trialdat;

rhits = params.trialid{2};
lhits = params.trialid{3};

trials = [rhits;lhits]; % using right and left hit trials only 

lfads = lfads(:,:,trials);

gauss = gauss(:,:,trials);

for i = 1:size(gauss,2) % for each clu
    gauss(:,i,:) = mySmooth(squeeze(gauss(:,i,:)),sm,sigma); % causal gaussian filter
end



lfads_recon = nan(numel(trials),1);
gauss_recon = nan(numel(trials),1);
for i = 1:numel(trials)
    test = i;
    train = trials~=trials(i);
    
    lfads_train = lfads(:,:,train);
    lfads_test = lfads(:,:,test);
    
    gauss_train = gauss(:,:,train);
    gauss_test = gauss(:,:,test);
    
    psth = mean(gauss_train,3);
    
    lfads_recon(i) = corr(lfads_test(:),psth(:));
    
    gauss_recon(i) = corr(gauss_test(:),psth(:));
    
end



%%

f = figure;
cols = [130, 176, 120 ; 106, 115, 171] ./ 255;
vs = violinplot([lfads_recon, gauss_recon],{'LFADS','Gaussian'},...
                'ViolinColor',cols,'EdgeColor',[1 1 1]);
ylabel('Neural Reconstruction, CC')
ylim([0,1])
ax = gca;
ax.FontSize = 20;
% 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/reconstruction';
% fn = 'recon_JEB7_2021-04-29_run12_gausssm_31_rhit_lhit';
% mysavefig(f,pth,fn)





%% visualize lfads factors
clear temp

lfads = dat.factors;

rhits = params.trialid{2};
lhits = params.trialid{3};

trials = [rhits;lhits]; % using right and left hit trials only 

temp{1} = lfads(:,:,rhits);
temp{2} = lfads(:,:,lhits);
lfads = temp; clear temp

%% trial-avg

clrs = getColors();

sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

f = figure;
alph = 0.5;
for i = 1:size(lfads{1},2)
    ax = nexttile; hold on;

    % lfads
    temp = squeeze(lfads{1}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(lfads{1},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs.rhit,'LineWidth',3},alph, ax);
    temp = squeeze(lfads{2}(:,i,:));
    stderr = std(temp,[],2) ./ sqrt(size(lfads{2},3));
    shadedErrorBar(obj.time, mean(temp,2), stderr, {'Color',clrs.lhit,'LineWidth',3},alph, ax);

    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    title(['Factor ' num2str(i)])
%     legend({'LFADS','Smoothing'},'Location','best')
    ax = gca;
    ax.FontSize = 20;
    xlim([obj.time(10) obj.time(end)])
    
    
end

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/lfads_factors';
% fn = [meta.anm '_' meta.date 'allfactors'  '_rlhits_lfads_' params.lfads_run ];
% mysavefig(f,pth,fn)

%% single trials


clrs = getColors();

sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

f = figure;
alph = 0.5;
for i = 1:size(lfads{1},2)
    ax = nexttile; hold on;

    % lfads
    temp = squeeze(lfads{1}(:,i,:));
    plot(obj.time,temp,'Color',clrs.rhit);
    temp = squeeze(lfads{2}(:,i,:));
    plot(obj.time,temp,'Color',clrs.lhit);

    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    title(['Factor ' num2str(i)])
%     legend({'LFADS','Smoothing'},'Location','best')
    ax = gca;
    ax.FontSize = 20;
    xlim([obj.time(10) obj.time(end)])
    
    
end


% 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/lfads_factors';
% fn = [meta.anm '_' meta.date 'allfactors_singletrials'  '_rlhits_lfads_' params.lfads_run ];
% mysavefig(f,pth,fn)




%% show condition-averaged psth, lfads single trials, gaussian smoothed single trials

clear temp

lfads = dat.rates;

gauss = obj.trialdat;

rhits = params.trialid{2};
lhits = params.trialid{3};


temp{1} = lfads(:,:,rhits);
temp{2} = lfads(:,:,lhits);


lfads = temp; clear temp

temp{1} = gauss(:,:,rhits);
temp{2} = gauss(:,:,lhits);


gauss = temp; clear temp



sm = 51; sigma = 51; winsize = sigma;
for j = size(gauss)
    for i = 1:size(gauss{j},2) % for each clu
        gauss{j}(:,i,:) = mySmooth(squeeze(gauss{j}(:,i,:)),sm,sigma); % causal gaussian filter
    end
end
% save same cells as above

close all
sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

clrs = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840],...
        [184, 146, 79]./255,[78, 186, 96]./255};

alph = 0.5;
% clus =  randsample(params.cluid,1);
clus = [21,34,54];
for i = 1:numel(clus)
    f = figure;
    f.Position= [-1666         165        1144         273];
    clf(f); 
    
    ix = find(params.cluid==clus(i));
    
    % gaussian psth
    ax = subplot(1,3,1); hold on
    for j = 1:numel(gauss)
        temp = squeeze(gauss{j}(:,ix,:));
        plot(obj.time,mean(temp,2),'Color',clrs{j},'LineWidth',4)
    end
    hold off
    xlim([obj.time(10) obj.time(end)])
    title('PSTH')
    xlabel('Time (s) from go cue')
    ax.YTick = [];
    ax = gca;
    ax.FontSize = 20;
    

    
    % lfads single trials
    ax = subplot(1,3,2); hold on
    for j = 1:numel(lfads)
        temp = squeeze(lfads{j}(:,ix,:));
        plot(obj.time,temp,'Color',clrs{j})
    end
    xlim([obj.time(10) obj.time(end)])
    title('LFADS')
    xlabel('Time (s) from go cue')
    ax.YTick = [];
    ax = gca;
    ax.FontSize = 20;
    hold off
    
    
    % gauss single trials
    ax = subplot(1,3,3); hold on
    for j = 1:numel(gauss)
        temp = squeeze(gauss{j}(:,ix,:));
        plot(obj.time,temp,'Color',clrs{j})
    end
    title(['Gaussian | sigma=' num2str(sigma) ' ms'])
    xlim([obj.time(10) obj.time(end)])
    xlabel('Time (s) from go cue')
    ax.YTick = [];
    ax = gca;
    ax.FontSize = 20;
    hold off
    
    
    pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/psth_lfads_gauss_compare';
    fn = [meta.anm '_' meta.date '_' params.lfads_run '_cell_' num2str(clus(i)) ];
    mysavefig(f,pth,fn)
    
    pause(10)
    
end












