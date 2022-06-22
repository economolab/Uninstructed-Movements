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
% meta = loadJEB4_ALMVideo(meta,datapth); % done
% meta = loadJEB5_ALMVideo(meta,datapth); % done
meta = loadJEB6_ALMVideo(meta,datapth); % done
meta = loadJEB7_ALMVideo(meta,datapth); % done
meta = loadEKH1_ALMVideo(meta,datapth); % done
meta = loadEKH3_ALMVideo(meta,datapth); % done
meta = loadJGR2_ALMVideo(meta,datapth); % done
meta = loadJGR3_ALMVideo(meta,datapth); % done


params.probe = [meta.probe]; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

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

%% Motion Energy

use = true(size(meta));
for i = 1:numel(meta)
    me(i) = loadMotionEnergy(objs{i},meta(i),params,1:objs{i}.bp.Ntrials); 
    if ~me(i).use
        use(i) = false;
    end
    disp('DONE');
    disp(' ');
end

params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);
meta = meta(use);
objs = objs(use);
me = me(use);

%%

clearvars -except meta params objs me


%% elsayed method single trials

clear rez

for i = 1:numel(meta)
    input_data = objs{i}.trialdat; % dat(i).factors   dat(i).rates (time,neurons,trials)
    rez(i) = elsayedNullandPotentSpace(objs{i},input_data,me(i),params);
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
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace';
% fn = 've_total';
% mysavefig(f,pth,fn);

violincols = [50, 168, 82; 168, 50, 142; 168, 50, 142; 50, 168, 82] ./ 225;
varexp_sep = [null_prep; null_move ; potent_move; potent_prep]';
f = figure; ax = axes(f);
vs = violinplot(varexp_sep,{'Null, Non-move','Null, Move', 'Potent,Move','Potent, Non-move'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1},  'ViolinColor', violincols);
ylabel('Normalized Variance Explained')
ylim([0,1])
ax = gca;
ax.FontSize = 20;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace';
% fn = 've_allsep';
% mysavefig(f,pth,fn);

%% plot projections

close all

sav = 0;

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
            tempdat = squeeze(temp(:,params.trialid{i}{j+1},dimix));
            means = mean(tempdat,2);
            stderr = std(tempdat,[],2) ./ numel(params.trialid{i}{j+1});
            shadedErrorBar(objs{1}.time,means,stderr,{'Color',clrs{j},'LineWidth',lw},alph, ax)
%             plot(objs{1}.time,means,'Color',clrs{j},'LineWidth',1)
        end
        title(['Potent ' num2str(dimix)])
        xlim([objs{1}.time(15),objs{1}.time(end)])
        
        align = mode(objs{1}.bp.ev.(params.alignEvent));
        sample = mode(objs{1}.bp.ev.sample) - align;
        delay = mode(objs{1}.bp.ev.delay) - align;
        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)
        
        ax.FontSize = 20;
        hold off;
    end
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace/potent';
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
            tempdat = squeeze(temp(:,params.trialid{i}{j+1},dimix));
            means = mean(tempdat,2);
            stderr = std(tempdat,[],2) ./ numel(params.trialid{i}{j+1});
            shadedErrorBar(objs{1}.time,means,stderr,{'Color',clrs{j},'LineWidth',lw},alph, ax)
%             plot(objs{1}.time,means,'Color',clrs{j},'LineWidth',1)
        end
        title(['Null ' num2str(dimix)])
        xlim([objs{1}.time(15),objs{1}.time(end)])
        
        align = mode(objs{1}.bp.ev.(params.alignEvent));
        sample = mode(objs{1}.bp.ev.sample) - align;
        delay = mode(objs{1}.bp.ev.delay) - align;
        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)
        
        ax.FontSize = 20;
        hold off;
    end
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace/null';
        fn = [meta(i).anm '_' meta(i).date];
        mysavefig(f,pth,fn);
        pause(3)
    end
    
    
    
end


%% activity modes

clear cdrez times

for i = 1:numel(rez)
    [cdrez(i),times] = cdNullSpace_binnedRates(rez(i),objs{i},params,i);
end

rez = cdrez;

% modes.varexp = activityModes_varexp(modes,rez);

%% null space cds

% concatenate latents, find mean and stderror
cdEarly{1} = rez(1).cd.null.cdEarly_latent(:,params.trialid{1}{2});
cdEarly{2} = rez(1).cd.null.cdEarly_latent(:,params.trialid{1}{3});
cdLate{1} = rez(1).cd.null.cdLate_latent(:,params.trialid{1}{2});
cdLate{2} = rez(1).cd.null.cdLate_latent(:,params.trialid{1}{3});
cdGo{1} = rez(1).cd.null.cdGo_latent(:,params.trialid{1}{2});
cdGo{2} = rez(1).cd.null.cdGo_latent(:,params.trialid{1}{3});
for i = 2:numel(rez)
    for j = 1:2
        cdEarly{j} = cat(2,cdEarly{j},rez(i).cd.null.cdEarly_latent(:,params.trialid{i}{j+1}));
        cdLate{j} = cat(2,cdLate{j},rez(i).cd.null.cdLate_latent(:,params.trialid{i}{j+1}));
        cdGo{j} = cat(2,cdGo{j},rez(i).cd.null.cdGo_latent(:,params.trialid{i}{j+1}));
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

sample = mode(objs{1}.bp.ev.sample - objs{1}.bp.ev.(params.alignEvent));
delay = mode(objs{1}.bp.ev.delay - objs{1}.bp.ev.(params.alignEvent));

fns = patternMatchCellArray(fieldnames(rez(1).cd.null),{'mode'},'all');

sav = 0;
for i = 1:numel(fns)
    f(i) = figure; ax = axes(f(i)); hold on
    tempmean = mySmooth(eval([fns{i}(1:end-5) '_latent_mean']),sm);
    temperror = mySmooth(eval([fns{i}(1:end-5) '_latent_error']),sm);
    shadedErrorBar(objs{1}.time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(objs{1}.time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([objs{1}.time(10);objs{1}.time(end)])
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
    shadetimes = objs{1}.time(times.(timefns{ix}));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
%     ylim(ylims);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace_binnedRates/null';
        fn = [fns{i}];
        mysavefig(f(i),pth,fn);
    end
    
end

% clear selectivity
% 
% sumsqselectivity.total = zeros(numel(objs{1}.time),numel(rez));
% sumsqselectivity.early = zeros(numel(objs{1}.time),numel(rez));
% sumsqselectivity.late = zeros(numel(objs{1}.time),numel(rez));
% sumsqselectivity.go = zeros(numel(objs{1}.time),numel(rez));
% for i = 1:numel(rez)
%     temp = squeeze(nanmean(gpfa(i).gpfalatents(:,:,params(i).trialid{2}),3));
%     rez(i).gpfa(:,:,1) = temp;
%     temp = squeeze(nanmean(gpfa(i).gpfalatents(:,:,params(i).trialid{3}),3));
%     rez(i).gpfa(:,:,2) = temp;
%     
%     selectivity.total = rez(i).gpfa(:,:,1) - rez(i).gpfa(:,:,2);
%     selectivity.early = mean(rez(i).cd.null.cdEarly_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.null.cdEarly_latent(:,params(i).trialid{3}),2);
%     selectivity.late  = mean(rez(i).cd.null.cdLate_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.null.cdLate_latent(:,params(i).trialid{3}),2);
%     selectivity.go    = mean(rez(i).cd.null.cdGo_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.null.cdGo_latent(:,params(i).trialid{3}),2);
%     selfns = fieldnames(selectivity);
%     for j = 1:numel(selfns)
%         selectivity.(selfns{j})(isnan(selectivity.(selfns{j}))) = 0;
%         sumsqselectivity.(selfns{j})(:,i) = sum(selectivity.(selfns{j}).^2,2);
%     end
% end
% 
% nullselexp_early = mean(sumsqselectivity.total ./ (sumsqselectivity.early),1)';
% nullselexp_late = mean(sumsqselectivity.total ./ (sumsqselectivity.late),1)';
% nullselexp_go = mean(sumsqselectivity.total ./ (sumsqselectivity.go),1)';
% 
% nullselexp = sumsqselectivity.total ./ (sumsqselectivity.early+sumsqselectivity.late+sumsqselectivity.go);
% nullselexp_tot = mean(nullselexp,1)';



%% potent space cds

% concatenate latents, find mean and stderror

cdEarly{1} = rez(1).cd.potent.cdEarly_latent(:,params.trialid{1}{2});
cdEarly{2} = rez(1).cd.potent.cdEarly_latent(:,params.trialid{1}{3});
cdLate{1} = rez(1).cd.potent.cdLate_latent(:,params.trialid{1}{2});
cdLate{2} = rez(1).cd.potent.cdLate_latent(:,params.trialid{1}{3});
cdGo{1} = rez(1).cd.potent.cdGo_latent(:,params.trialid{1}{2});
cdGo{2} = rez(1).cd.potent.cdGo_latent(:,params.trialid{1}{3});
for i = 2:numel(rez)
    for j = 1:2
        cdEarly{j} = cat(2,cdEarly{j},rez(i).cd.potent.cdEarly_latent(:,params.trialid{i}{j+1}));
        cdLate{j} = cat(2,cdLate{j},rez(i).cd.potent.cdLate_latent(:,params.trialid{i}{j+1}));
        cdGo{j} = cat(2,cdGo{j},rez(i).cd.potent.cdGo_latent(:,params.trialid{i}{j+1}));
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

sample = mode(objs{1}.bp.ev.sample - objs{1}.bp.ev.(params.alignEvent));
delay = mode(objs{1}.bp.ev.delay - objs{1}.bp.ev.(params.alignEvent));

fns = patternMatchCellArray(fieldnames(rez(1).cd.potent),{'mode'},'all');

sav = 0;
for i = 1:numel(fns)
    f(i) = figure; ax = axes(f(i)); hold on
    tempmean = mySmooth(eval([fns{i}(1:end-5) '_latent_mean']),sm);
    temperror = mySmooth(eval([fns{i}(1:end-5) '_latent_error']),sm);
    shadedErrorBar(objs{1}.time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(objs{1}.time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([objs{1}.time(10);objs{1}.time(end)])
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
    shadetimes = objs{1}.time(times.(timefns{ix}));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
%     ylim(ylims);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace_binnedRates/potent';
        fn = [fns{i}];
        mysavefig(f(i),pth,fn);
    end
    
end

% clear selectivity
% 
% sumsqselectivity.total = zeros(numel(objs{1}.time),numel(rez));
% sumsqselectivity.early = zeros(numel(objs{1}.time),numel(rez));
% sumsqselectivity.late = zeros(numel(objs{1}.time),numel(rez));
% sumsqselectivity.go = zeros(numel(objs{1}.time),numel(rez));
% for i = 1:numel(rez)
%     temp = squeeze(nanmean(gpfa(i).gpfalatents(:,:,params(i).trialid{2}),3));
%     rez(i).gpfa(:,:,1) = temp;
%     temp = squeeze(nanmean(gpfa(i).gpfalatents(:,:,params(i).trialid{3}),3));
%     rez(i).gpfa(:,:,2) = temp;
%     
%     selectivity.total = rez(i).gpfa(:,:,1) - rez(i).gpfa(:,:,2);
%     selectivity.early = mean(rez(i).cd.potent.cdEarly_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.potent.cdEarly_latent(:,params(i).trialid{3}),2);
%     selectivity.late  = mean(rez(i).cd.potent.cdLate_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.potent.cdLate_latent(:,params(i).trialid{3}),2);
%     selectivity.go    = mean(rez(i).cd.potent.cdGo_latent(:,params(i).trialid{2}),2) - mean(rez(i).cd.potent.cdGo_latent(:,params(i).trialid{3}),2);
%     selfns = fieldnames(selectivity);
%     for j = 1:numel(selfns)
%         selectivity.(selfns{j})(isnan(selectivity.(selfns{j}))) = 0;
%         sumsqselectivity.(selfns{j})(:,i) = sum(selectivity.(selfns{j}).^2,2);
%     end
% end
% 
% potentselexp_early = mean(sumsqselectivity.total ./ (sumsqselectivity.early),1)';
% potentselexp_late = mean(sumsqselectivity.total ./ (sumsqselectivity.late),1)';
% potentselexp_go = mean(sumsqselectivity.total ./ (sumsqselectivity.go),1)';
% 
% potentselexp = sumsqselectivity.total ./ (sumsqselectivity.early+sumsqselectivity.late+sumsqselectivity.go);
% potentselexp_tot = mean(potentselexp,1)';

%%
% close all
% 
% sav = 0;
% 
% violincols = [50, 168, 82; 168, 50, 142] ./ 255;
% f = figure; ax = axes(f);
% vs = violinplot([1./nullselexp_tot 1./potentselexp_tot],{'Null','Potent'},...
%     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
% ylabel('Selectivity Ratio (Sum / Total)')
% % ylim([0,1])
% ax = gca;
% ax.FontSize = 25;
% if sav
%     pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace/';
%     fn = 'selexp_null_potent';
%     mysavefig(f,pth,fn);
% end
% 
% violincols = [50, 168, 82; 168, 50, 142] ./ 255;
% f = figure; ax = axes(f);
% vs = violinplot([1./nullselexp_early 1./potentselexp_early],{'Null','Potent'},...
%     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
% ylabel('Selectivity Ratio (Early / Total)')
% % ylim([0,1])
% ax = gca;
% ax.FontSize = 25;
% if sav
%     pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace/';
%     fn = 'selexp_early';
%     mysavefig(f,pth,fn);
% end
% 
% violincols = [50, 168, 82; 168, 50, 142] ./ 255;
% f = figure; ax = axes(f);
% vs = violinplot([1./nullselexp_late 1./potentselexp_late],{'Null','Potent'},...
%     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
% ylabel('Selectivity Ratio (Late / Total)')
% % ylim([0,1])
% ax = gca;
% ax.FontSize = 25;
% if sav
%     pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace/';
%     fn = 'selexp_late';
%     mysavefig(f,pth,fn);
% end
% 
% violincols = [50, 168, 82; 168, 50, 142] ./ 255;
% f = figure; ax = axes(f);
% vs = violinplot([1./nullselexp_go 1./potentselexp_go],{'Null','Potent'},...
%     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor', violincols);
% ylabel('Selectivity Ratio (Go / Total)')
% % ylim([0,1])
% ax = gca;
% ax.FontSize = 25;
% if sav
%     pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace/';
%     fn = 'selexp_go';
%     mysavefig(f,pth,fn);
% end



