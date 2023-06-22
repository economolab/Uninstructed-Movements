clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'fig2/')))
rmpath(genpath(fullfile(utilspth,'figx/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))
rmpath(genpath(fullfile(utilspth,'musall2019/')))

% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                                              % all trials       (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % left hits, aw    (9)

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect'; % reflect, zeropad, none

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'Excellent','Great','Good','Fair'};

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth); % selectivity in ME
meta = loadEKH1_ALMVideo(meta,datapth); % selectivity in ME
meta = loadEKH3_ALMVideo(meta,datapth); % selectivity in ME
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth); % selectivity in ME % go cue is at 2.3 instead of 2.5 like all other sessions??
% meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);


% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

for sessix = 1:numel(meta)
    disp(['Loading kinematics - Session ' num2str(sessix) '/' num2str(numel(meta))])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end

%% CODING DIMENSIONS

clearvars -except obj meta params sel_corr_mat datapth me kin

% % 2afc (early, late, go)
cond2use = [2 3]; % left hit, right hit
cond2proj = [2 3 4 5 6 7 8 9];
rez_2afc = getCodingDimensions_2afc(obj,params,cond2use,cond2proj, 1);

% % aw (context mode)
cond2use = [6 7]; % hit 2afc, hit aw
cond2proj = [2 3 4 5 6 7 8 9];
rez_aw = getCodingDimensions_aw(obj,params,cond2use,cond2proj);


allrez = concatRezAcrossSessions(rez_2afc,rez_aw);

%% CONTEXT MODE ON SINGLE TRIALS

% CHANGE this to CDchoice/ramp %%%%%%%%%%%%%%%%

for sessix = 1:numel(meta)
    Wcontext = rez_aw(sessix).cd_mode_orth;
    trialdat{sessix} = tensorprod(obj(sessix).trialdat,Wcontext,2,1); % (time,trials)
end


%% DECODING PARAMETERS

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)

par.pre=6; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

% these parameters above are important for decoding accuracy. for actual
% analyses (like a final analysis to be included in a publication),
% you should vary these parameters only if you have a validation
% set that you won't test on until these parameters are set. Otherwise,
% there's a high risk of overfitting

% cross val folds
par.nFolds = 4;

% data sets
par.train = 0.7; % fraction of trials 
par.test = 1 - par.train;

% feature to use to decode
par.feats = kin(1).featLeg;
par.feats = {'tongue','motion','nose','jaw'};
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});
% par.feats = {'tongue_angle','tongue_length','motion_energy'};

% trials
par.cond2use = [6 7]; % CHANGE this for each CD prediction %%%%%%%%%%%%%%%%

par.regularize = 1; % if 0, linear regression. if 1, ridge regression

par.lambdas = logspace(-3,3,20); % regularization parameters to search over. don't need to change unless you want to make it more fine-grained. 

%% DECODING

close all

cmapfn = 'ContextColormap.mat'; % CHANGE probably
temp = load(fullfile('C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3\utils',cmapfn));
contextcmap = flipud(temp.ContextColormap);

for sessix = 1:numel(meta) % 7 is the example session for CDcontext
    disp([num2str(sessix) ' / ' num2str(numel(meta))])

    % predict 'trialdat' from par.feats (kinematic features of interest)

    % trials
    par.trials.all = cell2mat(params(sessix).trialid(par.cond2use)');
    % get delayed-response and water-cued trials so that we can match them
    % up with the test trials later
    par.trials.dr = params(sessix).trialid{par.cond2use(1)}; % CHANGE maybe %%%%%%%%%%%%%%%%%
    par.trials.wc = params(sessix).trialid{par.cond2use(2)}; % CHANGE maybe %%%%%%%%%%%%%%%%%

    % partition train and test
    nTrials = numel(par.trials.all);
    nTrain = floor(nTrials*par.train);
    par.trials.train = randsample(par.trials.all,nTrain,false);
    par.trials.test = par.trials.all(~ismember(par.trials.all,par.trials.train));

    % input data
    par.featix = find(ismember(kin(sessix).featLeg,par.feats));

    X.train = kin(sessix).dat(:,par.trials.train,par.featix); % (time,trials,feats)
    X.size.train = size(X.train);
    X.train = reshape(X.train, size(X.train,1)*size(X.train,2),size(X.train,3));

    X.test = kin(sessix).dat(:,par.trials.test,par.featix); % (time,trials,feats)
    X.size.test = size(X.test);
    X.test = reshape(X.test, size(X.test,1)*size(X.test,2),size(X.test,3));

    % reshape train and test data to account for prediction bin size
    X.train = reshapePredictors(X.train,par);
    X.test = reshapePredictors(X.test,par);

    % flatten inputs
    % if you're using a model with recurrence, don't flatten
    X.train = reshape(X.train,size(X.train,1),size(X.train,2)*size(X.train,3));
    X.test = reshape(X.test,size(X.test,1),size(X.test,2)*size(X.test,3));

    % output data
    Y.train = trialdat{sessix}(:,par.trials.train); % (time,trials);
    Y.size.train = size(Y.train);
    Y.train = reshape(Y.train, size(Y.train,1)*size(Y.train,2),size(Y.train,3));

    Y.test = trialdat{sessix}(:,par.trials.test);
    Y.size.test = size(Y.test);
    Y.test = reshape(Y.test, size(Y.test,1)*size(Y.test,2),size(Y.test,3));

    % standardize data
    % standardize both train and test sets using train set statistics
    % can also standardize using specific time points (presample for example)
    X.mu = mean(X.train,1,'omitnan');
    X.sigma = std(X.train,[],1,'omitnan');
    X.train = (X.train - X.mu) ./ X.sigma;
    if ~par.test==0
        X.test = (X.test - X.mu) ./ X.sigma;
    end

    Y.mu = mean(Y.train,1,'omitnan');
    Y.sigma = std(Y.train,[],1,'omitnan');
    Y.train = (Y.train - Y.mu);% ./ Y.sigma;
    if ~par.test==0
        Y.test = (Y.test - Y.mu);% ./ Y.sigma;
    end

    % fill missing values in kinematics
    X.train = fillmissing(X.train,'constant',0);
    Y.train = fillmissing(Y.train,'nearest');
    X.test = fillmissing(X.test,'constant',0);
    Y.test = fillmissing(Y.test,'nearest');

    for ilambda = 1:numel(par.lambdas)
        lambda = par.lambdas(ilambda);
        cvmdl{ilambda} = fitrlinear(X.train,Y.train,'Learner','svm','KFold',par.nFolds,'Regularization','ridge','Lambda',lambda);
        loss(ilambda) = cvmdl{ilambda}.kfoldLoss;
        % cvpred{ilambda} = kfoldPredict(cvmdl{ilambda});
    end
    
    [~,bestmodel] = min(loss);
    mdl = cvmdl{bestmodel}; % we now have the best lambda, and they trained cvmodel with that lambda,
    % can predict test data using best cv mdl
    for i = 1:par.nFolds
        mdl_loss(i) = mdl.Trained{i}.loss(X.train,Y.train);
    end
    [~,bestmodel] = min(mdl_loss);
    testmdl = mdl.Trained{bestmodel}; % we now have the best lambda, and the trained model with that lambda,
    
    pred = testmdl.predict(X.test);

    % the R-squared should be for the data we plot in the scatterplots (for
    % cd context it's the r2 for the ITI. Going to just use fitlm() for this since it also plots the 
    % scatterplot for us with the best fit line).

    y = reshape(Y.test,Y.size.test(1),Y.size.test(2)); % original input data (centered)
    yhat = reshape(pred,Y.size.test(1),Y.size.test(2)); % prediction

    sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.(params(1).alignEvent));
    delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.(params(1).alignEvent));
    gc = mode(obj(1).bp.ev.goCue) - mode(obj(1).bp.ev.(params(1).alignEvent));

    % find which trials were DR and which were WC in test set
    dr_trials = ismember(par.trials.test,par.trials.dr); %% CHANGE %%%%%%%%%%%%
    wc_trials = ismember(par.trials.test,par.trials.wc); %% CHANGE %%%%%%%%%%%%

    % neural data, split into trial type
    ydr = y(:,dr_trials); %% CHANGE %%%%%%%%%%%%
    ywc = y(:,wc_trials); %% CHANGE %%%%%%%%%%%%

    % predicted neural data, split into ground truth
    yhatdr = yhat(:,dr_trials); %% CHANGE %%%%%%%%%%%%
    yhatwc = yhat(:,wc_trials); %% CHANGE %%%%%%%%%%%%


    %% PLOTS

    % heatmaps


    cols = getColors;
    sm = 31;
    alph = 0.2;
    xlims = [-2.5 2.5];
    lw = 2;

    heatmapsm = 11;

    f = figure;
    f.Renderer = 'painters';
    ax = subplot(1,2,1);
    hold on;
    temp = mySmooth(cat(2,ydr,ywc),heatmapsm);

    imagesc(obj(sessix).time,1:size(y,2),temp');
    ax.YDir = "reverse";
    yline(size(ydr,2),'k-');
    xline(sample,'k--')
    xline(delay,'k--')
    xline(gc,'k--')
    colormap(contextcmap)
    title('data')
    c1 = colorbar;
    ylim(ax,[1 size(y,2)])
    % change tick direction to outside
    ax.TickDir = 'out';

    ax = subplot(1,2,2);

    temp = mySmooth(cat(2,yhatdr,yhatwc),heatmapsm);
    imagesc(obj(sessix).time,1:size(yhat,2),temp');
    ax.YDir = "reverse";
    yline(size(ydr,2),'k-','linewidth',0.1);
    xline(sample,'k--','linewidth',0.1)
    xline(delay,'k--','linewidth',0.1)
    xline(gc,'k--','linewidth',0.1)
    colormap(contextcmap)
    title('prediction')
    c2 = colorbar;
    ylim(ax,[1 size(y,2)])
    % change tick direction to outside
    ax.TickDir = 'out';


    %% trial-averaged plot


    mu.ydr = mean(ydr,2,'omitmissing');
    mu.ywc = mean(ywc,2,'omitmissing');
    mu.yhatdr = mean(yhatdr,2,'omitmissing');
    mu.yhatwc = mean(yhatwc,2,'omitmissing');

    sd.ydr = getCI(ydr); % 95% CI
    sd.ywc = getCI(ywc);
    sd.yhatdr = getCI(yhatdr);
    sd.yhatwc = getCI(yhatwc);


    fns = fieldnames(mu);
    for i = 1:numel(fns)
        mu.(fns{i}) = mySmooth(mu.(fns{i}),sm,'reflect');
        sd.(fns{i}) = mySmooth(sd.(fns{i}),sm,'reflect');
    end


    f = figure;
    f.Position = [644   483   338   231];
    ax = gca;
    f.Renderer = 'painters';
    ax = prettifyPlot(ax);
    hold on;
    shadedErrorBar(obj(1).time,mu.ydr,sd.ydr,{'Color',cols.afc,'LineWidth',lw,'LineStyle','-'},alph,ax)
    shadedErrorBar(obj(1).time,mu.ywc,sd.ywc,{'Color',cols.aw,'LineWidth',lw,'LineStyle','-'},alph,ax)
    shadedErrorBar(obj(1).time,mu.yhatdr,sd.yhatdr,{'Color',cols.afc,'LineWidth',lw,'LineStyle','--'},alph,ax)
    shadedErrorBar(obj(1).time,mu.yhatwc,sd.yhatwc,{'Color',cols.aw,'LineWidth',lw,'LineStyle','--'},alph,ax)
    % ax.FontSize = 10;
    xlim(xlims)
    xline(sample,'k--')
    xline(delay,'k--')
    xline(gc,'k--')

    ylabel('CD Context projection')
    xlabel('Time from go cue (s)')
    % title(t,['R^2 = ' num2str(round(r2,2))])  
    
    %% scatter plot of average iti activity

    % CHANGE - time points to use
    tix = [-2.4,-2.2]; % presample 
    % tix = [-2.4,2.4];
    ix = findTimeIX(obj(1).time,tix);
    ix = ix(1):ix(2);

    % mean across time points
    sdr = mean(ydr(ix,:),1,'omitmissing');
    shatdr = mean(yhatdr(ix,:),1,'omitmissing');
    swc = mean(ywc(ix,:),1,'omitmissing');
    shatwc = mean(yhatwc(ix,:),1,'omitmissing');

    % fit a linear model, use coefficients to plot line and get var exp
    xxx = cat(2,sdr,swc);
    yyy = cat(2,shatdr,shatwc);
    mdl_ = fitlm(xxx,yyy);
    r2(sessix) = mdl_.Rsquared.Ordinary;

    W = mdl_.Coefficients.Estimate; % intercept, x1; slope, x2


    f = figure;
    f.Position = [644   483   338   231];
    ax = gca;
    f.Renderer = 'painters';
    ax = prettifyPlot(ax);
    hold on;
    scatter(shatdr,sdr,30,'filled','MarkerEdgeColor','none','MarkerFaceColor',cols.afc);
    scatter(shatwc,swc,30,'filled','MarkerEdgeColor','none','MarkerFaceColor',cols.aw);
    xlabel('CD Context, prediction')
    ylabel('CD Context, data')
    thisr2 = round(r2(sessix),2);
    title(['R^2 = ' num2str(thisr2) ])
    xs = ax.XLim;
    xs = linspace(xs(1),xs(2),100); % plot best fit line
    rr =  xs.* W(2) + W(1);
    plot(xs,rr,'k--')
    axis square

end

%% r2 bar plot

f = figure;
f.Position = [748   394   193   272];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

dat = r2;


xs = 0;
for i = 1:numel(xs)
    b(i) = bar(xs(i),mean(dat,'omitmissing'));
    b(i).FaceColor = [0.5 0.5 0.5]; %[143, 49, 139] ./ 255; % 143, 49, 139
    b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 0.5;
    b(i).BarWidth = 0.7;
    vs(i) = scatter(xs(i)*ones(numel(dat),1),dat,50,'MarkerFaceColor','none',...
        'MarkerEdgeColor','k','XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,nanmean(dat),nanstd(dat),'LineStyle','none','Color','k','LineWidth',1)
end
ax.XTick = [];
ylabel("R^2")
ax = gca;
% ax.FontSize = 12;
























