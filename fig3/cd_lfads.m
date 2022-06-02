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


%% FACTORS , AVERAGED BY CONDITION

rates(:,:,1) = squeeze(mean(dat.rates(:,:,params.trialid{2}),3));
rates(:,:,2) = squeeze(mean(dat.rates(:,:,params.trialid{3}),3));

factors(:,:,1) = squeeze(mean(dat.factors(:,:,params.trialid{2}),3));
factors(:,:,2) = squeeze(mean(dat.factors(:,:,params.trialid{3}),3));


%% coding dimensions with factors

rez.time = obj.time;
rez.psth = factors;
rez.condition = params.condition;
rez.alignEvent = params.alignEvent;
rez.ev = obj.bp.ev;

%% cd early mode

% e1 = mode(rez.ev.delay) - 0.5 - mode(rez.ev.(params.alignEvent));
% e2 = mode(rez.ev.delay) - 0.1 - mode(rez.ev.(params.alignEvent));

e1 = mode(rez.ev.sample) + 0.4 - mode(rez.ev.(params.alignEvent));
e2 = mode(rez.ev.sample) + 0.8 - mode(rez.ev.(params.alignEvent));

times.early = rez.time>e1 & rez.time<e2;
tempdat = rez.psth(:,:,[1,2]);
mu = squeeze(mean(tempdat(times.early,:,:),1));
sd = squeeze(std(tempdat(times.early,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cdEarly_mode = cd;

%% cd late mode

e1 = mode(rez.ev.goCue) - 0.5 - mode(rez.ev.(params.alignEvent));
e2 = mode(rez.ev.goCue) - 0.1 - mode(rez.ev.(params.alignEvent));

times.late = rez.time>e1 & rez.time<e2;
tempdat = rez.psth(:,:,[1,2]);
mu = squeeze(mean(tempdat(times.late,:,:),1));
sd = squeeze(std(tempdat(times.late,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cdLate_mode = cd;


%% cd go mode

e1 = mode(rez.ev.goCue) + 0.02 - mode(rez.ev.(params.alignEvent));
e2 = mode(rez.ev.goCue) + 0.42 - mode(rez.ev.(params.alignEvent));

times.go = rez.time>e1 & rez.time<e2;
tempdat = rez.psth(:,:,[1,2]);
mu = squeeze(mean(tempdat(times.go,:,:),1));
sd = squeeze(std(tempdat(times.go,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cdGo_mode = cd;

%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(size(rez.psth,2),numel(fns));
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

% normTimes{1} = rez.time>e1 & rez.time<e2; % sample
% normTimes{2} = rez.time>-0.6 & rez.time<0; % delay
% normTimes{3} = rez.time>0 & rez.time<0.4; % go

cond = [1 2];
for i = 1:numel(fns)
    tempmode = rez.(fns{i});
    for j = 1:numel(cond)
        c = cond(j);
        
        tempdat = rez.psth(:,:,c)*rez.(fns{i});
        
%         normfactor = abs(nanmean(tempdat(normTimes{i})));
        normfactor = 1;
        
        rez.([fns{i}(1:end-5) '_latent'])(:,j) = tempdat ./ normfactor;
    end
end

clear cond

%% variance explained

for i = 1:numel(fns)
    psth = rez.psth;
    datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
    datacov(isnan(datacov)) = 0;
    eigsum = sum(eig(datacov));
    rez.varexp.(fns{i}(1:end-5)) = var_proj(rez.(fns{i}), datacov, eigsum);
end


%% concatenate latents, find mean and stderror

cdEarly{1} = rez(1).cdEarly_latent(:,1);
cdEarly{2} = rez(1).cdEarly_latent(:,2);
cdLate{1} = rez(1).cdLate_latent(:,1);
cdLate{2} = rez(1).cdLate_latent(:,2);
cdGo{1} = rez(1).cdGo_latent(:,1);
cdGo{2} = rez(1).cdGo_latent(:,2);
for i = 2:numel(rez)
    for j = 1:2
        cdEarly{j} = cat(2,cdEarly{j},rez(i).cdEarly_latent(:,j));
        cdLate{j} = cat(2,cdLate{j},rez(i).cdLate_latent(:,j));
        cdGo{j} = cat(2,cdGo{j},rez(i).cdGo_latent(:,j));
    end
end

cdEarly_latent_mean = [nanmean(cdEarly{1},2) nanmean(cdEarly{2},2)];
cdLate_latent_mean = [nanmean(cdLate{1},2) nanmean(cdLate{2},2)];
cdGo_latent_mean = [nanmean(cdGo{1},2) nanmean(cdGo{2},2)];

cdEarly_latent_error = [nanstd(cdEarly{1},[],2) nanstd(cdEarly{2},[],2)] ./ numel(rez); % std error
cdLate_latent_error = [nanstd(cdLate{1},[],2) nanstd(cdLate{2},[],2)] ./ numel(rez); % std error
cdGo_latent_error = [nanstd(cdGo{1},[],2) nanstd(cdGo{2},[],2)] ./ numel(rez); % std error

%%
close all
clrs = getColors();
lw = 6;
alph = 0.5;

sample = mode(rez(1).ev.sample - rez(1).ev.(params(1).alignEvent));
delay = mode(rez(1).ev.delay - rez(1).ev.(params(1).alignEvent));

sav = 0;
for i = 1:numel(fns)
    f(i) = figure; ax = axes(f(i)); hold on
    tempmean = eval([fns{i}(1:end-5) '_latent_mean']);
    temperror = eval([fns{i}(1:end-5) '_latent_error']);
    shadedErrorBar(rez(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(rez(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([rez(1).time(15);rez(1).time(end)])
    ylims = [min(min(tempmean))-1, max(max(tempmean))+1];
    ylim(ylims);
    
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
    
    ylim(ylims);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/cd';
        fn = [fns{i} '_anmList1_sessionList1_w_excludedsessions_sm_' num2str(params.smooth)];
        mysavefig(f(i),pth,fn);
    end
    
end


%% selectivity

sumsqselectivity.total = zeros(numel(rez(1).time),numel(rez));
sumsqselectivity.early = zeros(numel(rez(1).time),numel(rez));
sumsqselectivity.late = zeros(numel(rez(1).time),numel(rez));
sumsqselectivity.go = zeros(numel(rez(1).time),numel(rez));
for i = 1:numel(rez)
    selectivity.total = rez(i).psth(:,:,1) - rez(i).psth(:,:,2);
    selectivity.early = rez(i).cdEarly_latent(:,1) - rez(i).cdEarly_latent(:,2);
    selectivity.late  = rez(i).cdLate_latent(:,1) - rez(i).cdLate_latent(:,2);
    selectivity.go    = rez(i).cdGo_latent(:,1) - rez(i).cdGo_latent(:,2);
    selfns = fieldnames(selectivity);
    for j = 1:numel(selfns)
        selectivity.(selfns{j})(isnan(selectivity.(selfns{j}))) = 0;
        sumsqselectivity.(selfns{j})(:,i) = sum(selectivity.(selfns{j}).^2,2);
    end
end

% selectivity explained
% calculate variance of selectivity matrix within each of the modes
clear selectivity
% only calculate variance explained between start of cdEarly and end of
% cdGo- all other selectivity is variance that we are not attempting to
% explain with the three CDs
ix1 = find(times.early,1,'first');
ix2 = find(times.go,1,'last');
for i = 1:numel(rez)
    selectivity(i).total = rez(i).psth(ix1:ix2,:,1) - rez(i).psth(ix1:ix2,:,2);
    selectivity(i).early = rez(i).cdEarly_latent(ix1:ix2,1) - rez(i).cdEarly_latent(ix1:ix2,2);
    selectivity(i).late  = rez(i).cdLate_latent(ix1:ix2,1) - rez(i).cdLate_latent(ix1:ix2,2);
    selectivity(i).go    = rez(i).cdGo_latent(ix1:ix2,1) - rez(i).cdGo_latent(ix1:ix2,2);
    
    varianceInTotalSelectivity = trace(cov(selectivity(i).total));
    selexp.early(i) = trace(cov(selectivity(i).early)) ./ varianceInTotalSelectivity;
    selexp.late(i) = trace(cov(selectivity(i).late)) ./ varianceInTotalSelectivity;
    selexp.go(i) = trace(cov(selectivity(i).go)) ./ varianceInTotalSelectivity;
    selexp.sum(i) = selexp.early(i) + selexp.late(i)  + selexp.go(i);
end
selectivityExplained = [selexp.early' selexp.late' selexp.go' selexp.sum'];

%%
close all
sample = mode(rez(1).ev.sample) - mode(rez(1).ev.(params(1).alignEvent));
delay  = mode(rez(1).ev.delay) - mode(rez(1).ev.(params(1).alignEvent));

clrs = {'k','b','g','m','c'};

lw = 4;
alph = 0.5;
f = figure; ax = axes(f); hold on;
summean = zeros(numel(rez(1).time),1);
sumstd = zeros(numel(rez(1).time),1);
for i = 1:numel(selfns)
    temp = sumsqselectivity.(selfns{i});
    tempmean = nanmean(temp,2);
    tempstd = nanstd(temp,[],2);
    shadedErrorBar(rez(1).time,tempmean,...
                   tempstd./(numel(rez)),...
                   {'Color',clrs{i},'LineWidth',lw},alph,ax);
   if ~strcmpi(selfns{i},'total')
       summean = summean + tempmean;
       sumstd = sumstd + tempstd;
   end
end
shadedErrorBar(rez(1).time,summean,...
                   sumstd./(numel(rez)),...
                   {'Color',clrs{end},'LineWidth',lw},alph,ax);

xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

xlabel('Time (s) from go cue')
ylabel('Squared Sum Selectivity')
% legend('Total selectivity','early','late','go','early + late + go')
xlim([rez(1).time(1)+0.2,rez(1).time(end)])
ax = gca;
ax.FontSize = 20;
 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/selectivity';
% fn = 'sqsumselectivity_hits_anmList1_sessionList1_psthsm_51_dt_0_005_lowFR_1';
% mysavefig(f,pth,fn);

% variance in selectivity explained
varfns = fieldnames(selexp);
f = figure; ax = axes(f);
vs = violinplot(selectivityExplained,{'early','late','go','sum'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1});
ylabel('Selectivity Explained')
ylim([0,1])
ax = gca;
ax.FontSize = 25;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/selectivity';
% fn = 'selexp_hits_anmList1_sessionList1_psthsm_51_dt_0_005_lowFR_1';
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



































