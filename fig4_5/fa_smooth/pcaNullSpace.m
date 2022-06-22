clear,clc,close all


addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta(1).anm = 'JEB7'; 
meta(1).date = '2021-04-29'; 

meta(1).anm = 'JEB7'; 
meta(1).date = '2021-04-30'; 

meta(end+1).anm = 'JEB6'; 
meta(end).date = '2021-04-18'; 

meta(end+1).anm = 'JGR2'; 
meta(end).date = '2021-11-16';

meta(end+1).anm = 'JGR2'; 
meta(end).date = '2021-11-17';

meta(end+1).anm = 'JGR3'; 
meta(end).date = '2021-11-18';

meta(end+1).anm = 'EKH3'; 
meta(end).date = '2021-08-11';

meta(end+1).anm = 'EKH3'; 
meta(end).date = '2021-08-07';

meta(end+1).anm = 'EKH1'; 
meta(end).date = '2021-08-07';




meta = assignDataPath(meta);

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
dfparams = getDefaultParams();
% params.probe = 2; % change default params.probe from '1' to '2'


% --SPECIFY METHOD OF DENOISING SINGLE TRIALS
% we need denoised, smooth single trial neural activity. We have two
% options:
% 1) load data that's ALREADY been passed through an lfads model
% 2) perform Factor Analysis on binned single trial data, followed by
%    smoothing
dfparams.lfads_or_fa = 'lfads'; % 'lfads' or 'fa'
dfparams.lfads_run = 'run1'; % 'run3' , leave empty to use most recent run
dfparams.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter
dfparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
dfparams.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is
assert(strcmpi(dfparams.full_or_reduced,'reduced'),'method to use full dimensional data doesnt exist, params.full_or_reduced should be set to `reduced`')

% --SPECIFY TIME POINTS AND LAG TO USE
dfparams.prep = [-2.5 -0.05]; % initial and final time points (seconds) defining prep epoch, relative to alignevent
dfparams.move = [-2.5 1.5];   % initial and final time points (seconds) defining move epoch, relative to alignevent
dfparams.advance_movement = 0.025; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

use = true(size(meta));
for i = 1:numel(meta)
    disp(['Loading data for ' meta(i).anm ' ' meta(i).date]);
%     [meta(i),params(i),obj(i),dat(i)] = getNeuralActivity(meta(i),dfparams);
    fa(i) = getFAData(meta(i),'run2');
    params(i) = fa(i).params;
    if isfield(fa(i).obj,'meta')
        fa(i).obj = rmfield(fa(i).obj,'meta');
    end
    obj(i) = fa(i).obj;
    me(i) = loadMotionEnergy(obj(i),meta(i),params(i),1:obj(i).bp.Ntrials); 
    if ~me(i).use
        use(i) = false;
    end
    disp('DONE');
    disp(' ');
end

fa = fa(use);
params = params(use);
meta = meta(use);
obj = obj(use);
me = me(use);

%%

clearvars -except meta params obj dat fa dfparams me


%% 

% We now have a struct, dat, that contains:
% - trials:        trial numbers in use
% - factors:       neural data that's been reduced using lfads or factor analysis (time,factors,trials)
% - rates:         denoised single trial neural data from lfads or smoothing (time,cells,trials)




%% pca method single trials


clear rez

for i = 1:numel(meta)
    input_data = fa(i).falatents; % dat(i).factors   dat(i).rates
    rez(i) = pcaNullandPotentSpace(obj(i),input_data,me(i),params(i));
end

%% variance explained plots
close all

null_total = zeros(size(rez));
potent_total = zeros(size(rez));
null_prep = zeros(size(rez));
null_move = zeros(size(rez));
potent_move = zeros(size(rez));
potent_prep = zeros(size(rez));
for i = 1:numel(rez)
    null_total(i) = rez(i).ve.null_total;
    potent_total(i) = rez(i).ve.potent_total;
    null_prep(i) = rez(i).ve.null_prep;
    null_move(i) = rez(i).ve.null_move;
    potent_move(i) = rez(i).ve.potent_move;
    potent_prep(i) = rez(i).ve.potent_prep;
end

violincols = [50, 168, 82; 168, 50, 142] ./ 225;
varexp_full = [null_total ; potent_total]';
f = figure; ax = axes(f);
vs = violinplot(varexp_full,{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols);
ylabel('Normalized Variance Explained (Whole Trial)')
ylim([0,1])
ax = gca;
ax.FontSize = 20;
% 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace';
% fn = 've_total';
% mysavefig(f,pth,fn);

violincols = [50, 168, 82; 168, 50, 142; 168, 50, 142; 50, 168, 82] ./ 225;
varexp_sep = [null_prep; null_move ; potent_move; potent_prep]';
varexp_sep(:,1) = varexp_sep(:,1) - 0.001*rand(numel(meta),1);
f = figure; ax = axes(f);
vs = violinplot(varexp_sep,{'Null, Non-move','Null, Move', 'Potent,Move','Potent, Non-move'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1},  'ViolinColor', violincols);
ylabel('Normalized Variance Explained')
ylim([0,1])
ax = gca;
ax.FontSize = 20;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace';
% fn = 've_allsep';
% mysavefig(f,pth,fn);

%% plot projections

close all
clear cols clrs

sav = 1;

cols = getColors();
clrs{1} = cols.rhit;
clrs{2} = cols.lhit;
lw = 3;
alph = 0.5;
for i = 1:numel(rez)
    %     optimization_plots(rez,obj,dat,params); % old
    
    temp = rez(i).N_potent;
    
    f = figure;
    f.Position = [-1143         -27         357         848];
    for dimix = 1:size(temp,3)
        ax = subplot(size(temp,3),1,dimix); hold on
        for j = 1:2
            tempdat = squeeze(temp(:,params(i).trialid{j+1},dimix));
            means = mean(tempdat,2);
            stderr = std(tempdat,[],2) ./ numel(params(i).trialid{j+1});
            shadedErrorBar(obj(1).time,means,stderr,{'Color',clrs{j},'LineWidth',lw},alph, ax)
%             plot(obj(1).time,means,'Color',clrs{j},'LineWidth',1)
        end
        title(['Potent ' num2str(dimix)])
        xlim([obj(1).time(15),obj(1).time(end)])
        
        align = mode(obj(i).bp.ev.(params(i).alignEvent));
        sample = mode(obj(i).bp.ev.sample) - align;
        delay = mode(obj(i).bp.ev.delay) - align;
        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)
        
        ax.FontSize = 20;
        hold off;
    end
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/potent';
        fn = [meta(i).anm '_' meta(i).date];
        mysavefig(f,pth,fn);
        pause(3)
    end
    
    temp = rez(i).N_null;
    
    f = figure;
    f.Position = [-1501         187         357         631];
    for dimix = 1:size(temp,3)
        ax = subplot(size(temp,3),1,dimix); hold on
        for j = 1:2
            tempdat = squeeze(temp(:,params(i).trialid{j+1},dimix));
            means = mean(tempdat,2);
            stderr = std(tempdat,[],2) ./ numel(params(i).trialid{j+1});
            shadedErrorBar(obj(1).time,means,stderr,{'Color',clrs{j},'LineWidth',lw},alph, ax)
%             plot(obj(1).time,means,'Color',clrs{j},'LineWidth',1)
        end
        title(['Null ' num2str(dimix)])
        xlim([obj(1).time(15),obj(1).time(end)])
        
        align = mode(obj(i).bp.ev.(params(i).alignEvent));
        sample = mode(obj(i).bp.ev.sample) - align;
        delay = mode(obj(i).bp.ev.delay) - align;
        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)
        
        ax.FontSize = 20;
        hold off;
    end
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/null';
        fn = [meta(i).anm '_' meta(i).date];
        mysavefig(f,pth,fn);
        pause(3)
    end
    
    
end


%% activity modes

clear cdrez times

for i = 1:numel(rez)
    [cdrez(i),times] = cdNullSpace_elsayed(rez(i),obj(i),params(i));
end

rez = cdrez;

% modes.varexp = activityModes_varexp(modes,rez);

%% null space cds

% concatenate latents, find mean and stderror
cdEarly{1} = rez(1).cd.null.cdEarly_latent(:,params(1).trialid{2});
cdEarly{2} = rez(1).cd.null.cdEarly_latent(:,params(1).trialid{3});
cdLate{1} = rez(1).cd.null.cdLate_latent(:,params(1).trialid{2});
cdLate{2} = rez(1).cd.null.cdLate_latent(:,params(1).trialid{3});
cdGo{1} = rez(1).cd.null.cdGo_latent(:,params(1).trialid{2});
cdGo{2} = rez(1).cd.null.cdGo_latent(:,params(1).trialid{3});
for i = 2:numel(rez)
    for j = 1:2
        cdEarly{j} = cat(2,cdEarly{j},rez(i).cd.null.cdEarly_latent(:,params(i).trialid{j+1}));
        cdLate{j} = cat(2,cdLate{j},rez(i).cd.null.cdLate_latent(:,params(i).trialid{j+1}));
        cdGo{j} = cat(2,cdGo{j},rez(i).cd.null.cdGo_latent(:,params(i).trialid{j+1}));
    end
end

cdEarly_latent_mean = [nanmean(cdEarly{1},2) nanmean(cdEarly{2},2)];
cdLate_latent_mean = [nanmean(cdLate{1},2) nanmean(cdLate{2},2)];
cdGo_latent_mean = [nanmean(cdGo{1},2) nanmean(cdGo{2},2)];

cdEarly_latent_error = [nanstd(cdEarly{1},[],2) nanstd(cdEarly{2},[],2)] ./ (numel(rez)); % std error
cdLate_latent_error = [nanstd(cdLate{1},[],2) nanstd(cdLate{2},[],2)] ./ numel(rez); % std error
cdGo_latent_error = [nanstd(cdGo{1},[],2) nanstd(cdGo{2},[],2)] ./ numel(rez); % std error


close all
clrs = getColors();
lw = 6;
alph = 0.5;

sm = 1;

sample = mode(obj(1).bp.ev.sample - obj(1).bp.ev.(params(1).alignEvent));
delay = mode(obj(1).bp.ev.delay - obj(1).bp.ev.(params(1).alignEvent));

fns = patternMatchCellArray(fieldnames(rez(1).cd.null),{'mode'},'all');

sav = 1;

for i = 1:numel(fns)
    f(i) = figure; ax = axes(f(i)); hold on
    tempmean = mySmooth(eval([fns{i}(1:end-5) '_latent_mean']),sm);
    temperror = mySmooth(eval([fns{i}(1:end-5) '_latent_error']),sm);
    shadedErrorBar(obj(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(obj(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([obj(1).time(10);obj(1).time(end)])
%     ylims = [min(min(tempmean)), max(max(tempmean))+5];
%     ylim(ylims);
    
    title(fns{i},'Interpreter','none')
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax = gca;
    ax.FontSize = 40;
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    curmodename = fns{i};
    timefns = fieldnames(times);
    mask = strfind(timefns,lower(curmodename(3:end-5)));
    ix = cellfun(@(x) isempty(x),mask,'UniformOutput',false);
    ix = ~cell2mat(ix);
    shadetimes = obj(1).time(times.(timefns{ix}));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
%     ylim(ylims);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace_fa/null/';
        fn = [fns{i}];
        mysavefig(f(i),pth,fn);
    end
    
end

clear selectivity

sumsqselectivity.total = zeros(numel(obj(1).time),numel(rez));
sumsqselectivity.early = zeros(numel(obj(1).time),numel(rez));
sumsqselectivity.late = zeros(numel(obj(1).time),numel(rez));
sumsqselectivity.go = zeros(numel(obj(1).time),numel(rez));
for i = 1:numel(rez)
    temp = squeeze(nanmean(fa(i).falatents(:,:,params(i).trialid{2}),3));
    rez(i).fa(:,:,1) = temp;
    temp = squeeze(nanmean(fa(i).falatents(:,:,params(i).trialid{3}),3));
    rez(i).fa(:,:,2) = temp;
    
    selectivity.total = rez(i).fa(:,:,1) - rez(i).fa(:,:,2);
    selectivity.early = mean(rez(i).cd.null.cdEarly_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.null.cdEarly_latent(:,params(i).trialid{3}),2);
    selectivity.late  = mean(rez(i).cd.null.cdLate_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.null.cdLate_latent(:,params(i).trialid{3}),2);
    selectivity.go    = mean(rez(i).cd.null.cdGo_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.null.cdGo_latent(:,params(i).trialid{3}),2);
    selfns = fieldnames(selectivity);
    for j = 1:numel(selfns)
        selectivity.(selfns{j})(isnan(selectivity.(selfns{j}))) = 0;
        sumsqselectivity.(selfns{j})(:,i) = sum(selectivity.(selfns{j}).^2,2);
    end
end

nullselexp_early = mean(sumsqselectivity.total ./ (sumsqselectivity.early),1)';
nullselexp_late = mean(sumsqselectivity.total ./ (sumsqselectivity.late),1)';
nullselexp_go = mean(sumsqselectivity.total ./ (sumsqselectivity.go),1)';

nullselexp = sumsqselectivity.total ./ (sumsqselectivity.early+sumsqselectivity.late+sumsqselectivity.go);
nullselexp_tot = mean(nullselexp,1)';



%% potent space cds

% concatenate latents, find mean and stderror

cdEarly{1} = rez(1).cd.potent.cdEarly_latent(:,params(1).trialid{2});
cdEarly{2} = rez(1).cd.potent.cdEarly_latent(:,params(1).trialid{3});
cdLate{1} = rez(1).cd.potent.cdLate_latent(:,params(1).trialid{2});
cdLate{2} = rez(1).cd.potent.cdLate_latent(:,params(1).trialid{3});
cdGo{1} = rez(1).cd.potent.cdGo_latent(:,params(1).trialid{2});
cdGo{2} = rez(1).cd.potent.cdGo_latent(:,params(1).trialid{3});
for i = 2:numel(rez)
    for j = 1:2
        cdEarly{j} = cat(2,cdEarly{j},rez(i).cd.potent.cdEarly_latent(:,params(i).trialid{j+1}));
        cdLate{j} = cat(2,cdLate{j},rez(i).cd.potent.cdLate_latent(:,params(i).trialid{j+1}));
        cdGo{j} = cat(2,cdGo{j},rez(i).cd.potent.cdGo_latent(:,params(i).trialid{j+1}));
    end
end

cdEarly_latent_mean = [nanmean(cdEarly{1},2) nanmean(cdEarly{2},2)];
cdLate_latent_mean = [nanmean(cdLate{1},2) nanmean(cdLate{2},2)];
cdGo_latent_mean = [nanmean(cdGo{1},2) nanmean(cdGo{2},2)];

cdEarly_latent_error = [nanstd(cdEarly{1},[],2) nanstd(cdEarly{2},[],2)] ./ numel(rez); % std error
cdLate_latent_error = [nanstd(cdLate{1},[],2) nanstd(cdLate{2},[],2)] ./ numel(rez); % std error
cdGo_latent_error = [nanstd(cdGo{1},[],2) nanstd(cdGo{2},[],2)] ./ numel(rez); % std error


close all
clrs = getColors();
lw = 6;
alph = 0.5;

sm = 1;

sample = mode(obj(1).bp.ev.sample - obj(1).bp.ev.(params(1).alignEvent));
delay = mode(obj(1).bp.ev.delay - obj(1).bp.ev.(params(1).alignEvent));

fns = patternMatchCellArray(fieldnames(rez(1).cd.potent),{'mode'},'all');

sav = 1;


for i = 1:numel(fns)
    f(i) = figure; ax = axes(f(i)); hold on
    tempmean = mySmooth(eval([fns{i}(1:end-5) '_latent_mean']),sm);
    temperror = mySmooth(eval([fns{i}(1:end-5) '_latent_error']),sm);
    shadedErrorBar(obj(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(obj(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([obj(1).time(10);obj(1).time(end)])
%     ylims = [min(min(tempmean)), max(max(tempmean))+5];
%     ylim(ylims);
    
    title(fns{i},'Interpreter','none')
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax = gca;
    ax.FontSize = 40;
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    curmodename = fns{i};
    timefns = fieldnames(times);
    mask = strfind(timefns,lower(curmodename(3:end-5)));
    ix = cellfun(@(x) isempty(x),mask,'UniformOutput',false);
    ix = ~cell2mat(ix);
    shadetimes = obj(1).time(times.(timefns{ix}));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
    tempylim = ax.YLim;
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    ax.YLim = tempylim;
    
%     ylim(ylims);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace_fa/potent/';
        fn = [fns{i}];
        mysavefig(f(i),pth,fn);
    end
    
end

clear selectivity

sumsqselectivity.total = zeros(numel(obj(1).time),numel(rez));
sumsqselectivity.early = zeros(numel(obj(1).time),numel(rez));
sumsqselectivity.late = zeros(numel(obj(1).time),numel(rez));
sumsqselectivity.go = zeros(numel(obj(1).time),numel(rez));
for i = 1:numel(rez)
    temp = squeeze(nanmean(fa(i).falatents(:,:,params(i).trialid{2}),3));
    rez(i).fa(:,:,1) = temp;
    temp = squeeze(nanmean(fa(i).falatents(:,:,params(i).trialid{3}),3));
    rez(i).fa(:,:,2) = temp;
    
    selectivity.total = rez(i).fa(:,:,1) - rez(i).fa(:,:,2);
    selectivity.early = mean(rez(i).cd.potent.cdEarly_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.potent.cdEarly_latent(:,params(i).trialid{3}),2);
    selectivity.late  = mean(rez(i).cd.potent.cdLate_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.potent.cdLate_latent(:,params(i).trialid{3}),2);
    selectivity.go    = mean(rez(i).cd.potent.cdGo_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.potent.cdGo_latent(:,params(i).trialid{3}),2);
    selfns = fieldnames(selectivity);
    for j = 1:numel(selfns)
        selectivity.(selfns{j})(isnan(selectivity.(selfns{j}))) = 0;
        sumsqselectivity.(selfns{j})(:,i) = sum(selectivity.(selfns{j}).^2,2);
    end
end

potentselexp_early = mean(sumsqselectivity.total ./ (sumsqselectivity.early),1)';
potentselexp_late = mean(sumsqselectivity.total ./ (sumsqselectivity.late),1)';
potentselexp_go = mean(sumsqselectivity.total ./ (sumsqselectivity.go),1)';

potentselexp = sumsqselectivity.total ./ (sumsqselectivity.early+sumsqselectivity.late+sumsqselectivity.go);
potentselexp_tot = mean(potentselexp,1)';

%%
close all

sav = 1;

violincols = [50, 168, 82; 168, 50, 142] ./ 255;
f = figure; ax = axes(f);
vs = violinplot([1./nullselexp_tot 1./potentselexp_tot],{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
ylabel('Selectivity Ratio (Sum / Total)')
% ylim([0,1])
ax = gca;
ax.FontSize = 25;
if sav
    pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/cdve';
    fn = 'selexp_null_potent';
    mysavefig(f,pth,fn);
end

violincols = [50, 168, 82; 168, 50, 142] ./ 255;
f = figure; ax = axes(f);
vs = violinplot([1./nullselexp_early 1./potentselexp_early],{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
ylabel('Selectivity Ratio (Early / Total)')
% ylim([0,1])
ax = gca;
ax.FontSize = 25;
if sav
    pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/cdve';
    fn = 'selexp_early';
    mysavefig(f,pth,fn);
end

violincols = [50, 168, 82; 168, 50, 142] ./ 255;
f = figure; ax = axes(f);
vs = violinplot([1./nullselexp_late 1./potentselexp_late],{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
ylabel('Selectivity Ratio (Late / Total)')
% ylim([0,1])
ax = gca;
ax.FontSize = 25;
if sav
    pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/cdve';
    fn = 'selexp_late';
    mysavefig(f,pth,fn);
end

violincols = [50, 168, 82; 168, 50, 142] ./ 255;
f = figure; ax = axes(f);
vs = violinplot([1./nullselexp_go 1./potentselexp_go],{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
ylabel('Selectivity Ratio (Go / Total)')
% ylim([0,1])
ax = gca;
ax.FontSize = 25;
if sav
    pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/cdve';
    fn = 'selexp_go';
    mysavefig(f,pth,fn);
end


