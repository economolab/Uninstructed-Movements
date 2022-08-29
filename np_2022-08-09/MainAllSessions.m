clear,clc,close all
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 0.5; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off

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


params.advance_movement = 0.0;

%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
meta = loadJEB7_ALMVideo(meta,datapth); % done
meta = loadEKH3_ALMVideo(meta,datapth); % done
meta = loadJGR2_ALMVideo(meta,datapth); % done
meta = loadJGR3_ALMVideo(meta,datapth); % done
meta = loadJEB15_ALMVideo(meta,datapth);


params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

obj = loadObjs(meta);


for sessix = 1:numel(meta)
    for prbix = 1:numel(params.probe{sessix})
        disp('______________________________________________________')
        disp(['Processing data for session ' [meta(sessix).anm '_' meta(sessix).date ' | Probe' num2str(params.probe{sessix}(prbix))   ]])
        disp(' ')

        prbnum = params.probe{sessix}(prbix);

        [sessparams{sessix,prbix},sessobj{sessix,prbix}] = processData(obj(sessix),params,prbnum);
    end
end

% clean up sessparams and sessobj
for sessix = 1:numel(meta)
    params.trialid{sessix} = sessparams{sessix}.trialid;

    if numel(params.probe{sessix}) == 1
        params.cluid{sessix} = sessparams{sessix,1}.cluid{params.probe{sessix}};

        objs(sessix) = sessobj{sessix,1};
        objs(sessix).psth = objs(sessix).psth{params.probe{sessix}};
        objs(sessix).trialdat = objs(sessix).trialdat{params.probe{sessix}};
        objs(sessix).presampleFR = objs(sessix).presampleFR{params.probe{sessix}};
        objs(sessix).presampleSigma = objs(sessix).presampleSigma{params.probe{sessix}};
    elseif numel(params.probe{sessix}) == 2 % concatenate both probes worth of data

        params.cluid{sessix} = {sessparams{sessix,1}.cluid{params.probe{sessix}(1)}, sessparams{sessix,2}.cluid{params.probe{sessix}(2)} };

        objs(sessix) = sessobj{sessix,1};

        objs(sessix).psth = cat(2, objs(sessix).psth{1}, sessobj{sessix,2}.psth{2}); 
        objs(sessix).trialdat = cat(2, objs(sessix).trialdat{1}, sessobj{sessix,2}.trialdat{2}); 
        objs(sessix).presampleFR = cat(1, objs(sessix).presampleFR{1}, sessobj{sessix,2}.presampleFR{2}); 
        objs(sessix).presampleSigma = cat(1, objs(sessix).presampleSigma{1}, sessobj{sessix,2}.presampleSigma{2}); 
    end
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')


clear obj
obj = objs;
clear objs


%% convert params to a struct array

temp = params;
clear params
for sessix = 1:numel(meta)
    temp2 = rmfield(temp,{'probe','trialid','cluid'});
    temp2.probe = temp.probe{sessix};
    temp2.trialid = temp.trialid{sessix};
    temp2.cluid = temp.cluid{sessix};
    params(sessix) = temp2;
end

%% motion energy

for sessix = 1:numel(meta)
    fn = ['motionEnergy_' meta(sessix).anm '_' meta(sessix).date];
    fn = fullfile(datapth,'DataObjects',meta(sessix).anm,fn);
    temp = load(fn);
    temp = temp.me;
    if isstruct(temp.data)
        temp.data = temp.data.data;
        me(sessix) = temp;
    else
        me(sessix) = temp;
    end

    taxis = obj(sessix).time;
    alignTimes = obj(sessix).bp.ev.(params(sessix).alignEvent);
    tempme = zeros(numel(taxis),numel(me(sessix).data)); % (time,trials)
    for trix = 1:size(tempme,2)
        tempme(:,trix) = interp1(obj(sessix).traj{1}(trix).frameTimes-0.5-alignTimes(trix),me(sessix).data{trix},taxis);
    end
    me(sessix).data = tempme;

end

%% analysis type
np_type = 'optim'; % 'optim'  'kin'  'pca'  'optim_psth'   'optim_psth_me'

data_type = 'binned_rates';  % 'fa'   'binned_rates'   'gpfa'

%% null/potent

clearvars -except meta params obj fa gpfa dfparams me kin kinfeats kinfeats_reduced np_type data_type

switch np_type

    case 'optim'
        for i = 1:numel(meta)
            input_data = obj(i).trialdat;
            rez(i) = elsayedNullandPotentSpace(obj(i),input_data,me(i),params(i));
        end

    case 'kin'
        for i = 1:numel(meta)
            input_data = getInputData(obj(i),fa(i),gpfa(i),data_type,np_type,kinfeats_reduced{i});
            input_data.factors = input_data.neural;
            input_data.feats   = input_data.behav;
            rez(i) = kinematicNullAndPotentSpace(obj(i),input_data,dfparams,params(i),me(i));
        end

    case 'pca'
        for i = 1:numel(meta)
            input_data = getInputData(obj(i),fa(i),gpfa(i),data_type,np_type);
            input_data = input_data.neural;
            rez(i) = pcaNullandPotentSpace(obj(i),input_data,me(i),params(i));
        end
    
    case 'optim_psth'
        for i = 1:numel(meta)
            rez(i) = elsayedPSTH_NullandPotentSpace(obj(i),params(i),dfparams);
        end

    case 'optim_psth_me'
        for i = 1:numel(meta)
            rez(i) = elsayed_PSTH_ME_NullandPotentSpace(obj(i),me(i),params(i),dfparams,kin(i),kinfeats{i});
        end

end



% %% variance explained plots (TODO: NEED TO WORK ON THIS)
% 
% varexpPlots(rez);
% 
% %% plot projections
% 
% plotProjections(params,obj,rez)
% 


%% activity modes

clear cdrez times

if strcmpi(np_type,'optim_psth')

    for i = 1:numel(rez)
        [cdrez(i),times] = nullspace_cd_psth(rez(i),obj(i),params(i));
    end

    rez = cdrez;

    plotCD_psth(rez,obj,params,times,meta);

elseif ~strcmpi(np_type,'optim_psth_me')

    for i = 1:numel(rez)
        [cdrez(i),times] = nullspace_cd(rez(i),obj(i),params(i));
    end

    rez = cdrez;

    %     % plot mean CDs across sessions
    %     % null space cds
    %     nullSpaceCD(rez,obj,params,times)
    %     % potent space cds
    %     potentSpaceCD(rez,obj,params,times)

    % plot CDs for each session
    for i = 1:numel(rez)
%         plotSessionCD(rez(i),obj(i),params(i),times,meta(i))
%         % null space cds
%         plotSessionNullSpaceCD(rez(i),obj(i),params(i),times)
%         % potent space cds
%         plotSessionPotentSpaceCD(rez(i),obj(i),params(i),times)
    end

end


%% debugging functions

close all


% % plot single trial binned rates
% sessionix = 1;
% plotSingleTrialBinnedRates(obj,sessionix)

% % plot single trial neural factors (factor analysis)
% sessionix = 1;
% plotSingleTrialNeuralFactors(fa,obj,sessionix)


% % plot kinematics
% sessionix = 1;
% plotKinematics(kin,kinfeats,obj,sessionix);


% % plot reduced kinematics
% sessionix = 1;
% plotReducedKinematics(kinfeats_reduced,obj,sessionix);

% % plot motion energy
% sessionix = 1;
% plotMotionEnergy(obj,me,sessionix,kinfeats,kin);


% % plot motion energy with a potent dimension
% sessix = 2;
% dim = 1;
% plotME_PotentDim(obj,me,rez,sessix,dim,data_type)
% set(gca,'XLim',[256.3261  264.4115])


% % same for a null dimension
% sessix = 2;
% dim = 1;
% plotME_NullDim(obj,me,rez,sessix,dim,data_type)
% set(gca,'XLim',[256.3261  264.4115])





%% find top 4 dims that explain most variance (n/p projections)
for sessix = 1:numel(rez)
    eigsumNull = sum(eig(rez(sessix).covNull));
    eigsumPotent = sum(eig(rez(sessix).covPotent));
    for dimix = 1:rez(sessix).dPrep
        rez(sessix).ve.null.frac(dimix) = var_proj(rez(sessix).Qnull(:,dimix),rez(sessix).covNull,eigsumNull);
    end
    for dimix = 1:rez(sessix).dMove
        rez(sessix).ve.potent.frac(dimix) = var_proj(rez(sessix).Qpotent(:,dimix),rez(sessix).covPotent,eigsumPotent);
    end
    [~,temp] = sort(rez(sessix).ve.null.frac,'descend');
    rez(sessix).ve.null.top4 = temp(1:4);
    [~,temp] = sort(rez(sessix).ve.potent.frac,'descend');
    rez(sessix).ve.potent.top4 = temp(1:4);
end


%% prob of jaw movement onset of selectivity vs onset of selectivity in CDs

clrs = getColors();
lw = 1;
alph = 0.3;
sm = 21;

sav = 0;

%

for sessix = 1:numel(meta)

    motion{1} = me(sessix).data(:,params(sessix).trialid{2});
    motion{2} = me(sessix).data(:,params(sessix).trialid{3});


    % trial averaged motion energy
    f1 = figure; hold on;
    plot(obj(sessix).time,mean(motion{1},2),'Color',clrs.rhit,'LineWidth',2)
    plot(obj(sessix).time,mean(motion{2},2),'Color',clrs.lhit,'LineWidth',2)
    xlim([obj(sessix).time(10) obj(sessix).time(end)])
    ylim([ 5 max(max(mean(motion{1},2),mean(motion{2},2)))])
    xlabel('Time (s) from go cue')
    ylabel('Trial-avg motion energy')
    title([meta(sessix).anm '_' meta(sessix).date],'Interpreter','none')
    ax = gca;
    ax.FontSize = 13;
    hold off

    % single trial motion energy
    f2 = figure;
    temp = cat(2,motion{1},motion{2});
    imagesc(obj(sessix).time,1:size(temp,2),temp');
    colorbar; clim([0 80])
    xlabel('Time (s) from go cue')
    ylabel('Trial number')
    title([meta(sessix).anm '_' meta(sessix).date],'Interpreter','none')
    ax = gca;
    ax.FontSize = 13;
    hold off

    % single trial n/p projections, just first 4 dims
    f3 = figure;
    f3.Position = [680    56   992   922];
    for dimix = 1:4
        nulldim = rez(sessix).ve.null.top4(dimix);

        subplot(4,2,(dimix*2)-1)
        temp = cat(2,rez(sessix).N_null(:,params(sessix).trialid{2},nulldim), rez(sessix).N_null(:,params(sessix).trialid{3},nulldim));
        imagesc(obj(sessix).time, 1:size(temp,2), temp')
        if dimix == 1
            title('Null')
        elseif dimix==4
            xlabel('Time (s) from go cue')
            ax = gca;
            ax.FontSize = 13;
        end

        potentdim = rez(sessix).ve.potent.top4(dimix);
        subplot(4,2,(dimix*2))
        temp = cat(2,rez(sessix).N_potent(:,params(sessix).trialid{2},potentdim), rez(sessix).N_potent(:,params(sessix).trialid{3},potentdim));
        imagesc(obj(sessix).time, 1:size(temp,2), temp')
        if dimix == 1
            title('Potent')
        elseif dimix==4
            xlabel('Time (s) from go cue')
            ax = gca;
            ax.FontSize = 13;
        end
    end
    sgtitle([meta(sessix).anm '_' meta(sessix).date],'Interpreter','none')


    % single trial n/p coding direction projections
    npfn = fieldnames(rez(sessix).cd);
    cdfn = {'cdEarly_latent','cdLate_latent','cdGo_latent'};
    f4 = figure;
    f4.Position = [680    56   992   922];
    for cdix = 1:numel(cdfn)
        for npix = 1:numel(npfn)
            if npix == 1
                subplot(3,2,(cdix*2)-1)
            elseif npix == 2
                subplot(3,2,(cdix*2))
            end
            temp = cat(2,rez(sessix).cd.(npfn{npix}).(cdfn{cdix})(:,params(sessix).trialid{2}), rez(sessix).cd.(npfn{npix}).(cdfn{cdix})(:,params(sessix).trialid{3}));
            imagesc(obj(sessix).time, 1:size(temp,2), temp')


            if cdix == 1 && npix == 1
                title('CD Null')
            elseif cdix == 1 && npix == 2
                title('CD Potent')
            elseif cdix == 3 
                xlabel('Time (s) from go cue')
                ax = gca;
                ax.FontSize = 13;
            end
        end
    end
    sgtitle([meta(sessix).anm '_' meta(sessix).date],'Interpreter','none')

    % trial-averaged n/p coding direction projections
    npfn = fieldnames(rez(sessix).cd);
    cdfn = {'cdEarly_latent','cdLate_latent','cdGo_latent'};
    f5 = figure;
    f5.Position = [680    56   992   922];
    for cdix = 1:numel(cdfn)
        for npix = 1:numel(npfn)
            if npix == 1
                subplot(3,2,(cdix*2)-1)
            elseif npix == 2
                subplot(3,2,(cdix*2))
            end
            clear temp
            temp{1} = rez(sessix).cd.(npfn{npix}).(cdfn{cdix})(:,params(sessix).trialid{2});
            temp{2} = rez(sessix).cd.(npfn{npix}).(cdfn{cdix})(:,params(sessix).trialid{3});
            hold on;
            plot(obj(sessix).time, mean(temp{1},2),'Color',clrs.rhit,'LineWidth',2);
            plot(obj(sessix).time, mean(temp{2},2),'Color',clrs.lhit,'LineWidth',2);
            xlim([obj(sessix).time(10) obj(sessix).time(end)])

            if cdix == 1 && npix == 1
                title('CD Null')
            elseif cdix == 1 && npix == 2
                title('CD Potent')
            elseif cdix == 3 
                xlabel('Time (s) from go cue')
                ax = gca;
                ax.FontSize = 13;
            end
        end
    end
    sgtitle([meta(sessix).anm '_' meta(sessix).date],'Interpreter','none')

    if sav
        figpth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\np_2022-08-09\figs';
        savedir = [meta(sessix).anm '_' meta(sessix).date];
        pth = fullfile(figpth,savedir);
        fn = 'avgMotionEnergy';
        mysavefig(f1,pth,fn)
        pause(1)
        fn = 'motionEnergyHeatmap';
        mysavefig(f2,pth,fn)
        pause(1)
        fn = 'npHeatmap';
        mysavefig(f3,pth,fn)
        pause(1)
        fn = 'cdHeatmap';
        mysavefig(f4,pth,fn)
        pause(1)
        fn = 'avgCD';
        mysavefig(f5,pth,fn)
        pause(1)
    end


end



%% coding directions











































