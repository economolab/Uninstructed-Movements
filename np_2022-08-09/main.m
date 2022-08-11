clear,clc,close all
addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMALs/SESSIONs TO LOAD
% meta(1).anm = 'JEB7';
% meta(1).date = '2021-04-29';

% meta(1).anm = 'JEB7';
% meta(1).date = '2021-04-30';

meta(1).anm = 'JEB6';
meta(1).date = '2021-04-18';


meta = assignDataPath(meta); % CHANGE DATA PATH IN HERE

% ------------------------------------------------------------------ %
% get some parameters that are common to all sessions/analyses
dfparams = getDefaultParams();

% load data
[obj,fa,gpfa,params,me,kin,kinfeats,kinfeats_reduced] = loadProcessedData(meta);


%% analysis type
np_type = 'optim'; % 'optim'  'kin'  'pca'  'optim_psth'   'optim_psth_me'

data_type = 'binned_rates';  % 'fa'   'binned_rates'   'gpfa'

%% null/potent

clearvars -except meta params obj fa gpfa dfparams me kin kinfeats kinfeats_reduced np_type data_type

switch np_type

    case 'optim'
        for i = 1:numel(meta)
            input_data = getInputData(obj(i),fa(i),gpfa(i),data_type,np_type);
            input_data = input_data.neural;
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
        plotSessionCD(rez(i),obj(i),params(i),times,meta(i))
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





%% prob of jaw movement onset of selectivity vs onset of selectivity in CDs

% ix = find(ismember(kin.featLeg,'jaw_ydisp_view1'));
% ix = find(ismember(kin.featLeg,'jaw_yvel_view1'));
ix = find(ismember(kin.featLeg,'motion_energy'));

temp = kinfeats{1}(:,:,ix);

jaw{1} = temp(:,params.trialid{2});
jaw{2} = temp(:,params.trialid{3});

clrs = getColors();
lw = 1;
alph = 0.3;
sm = 21;

figure; hold on;
plot(obj.time,mean(jaw{1},2),'Color',clrs.rhit,'LineWidth',2)
plot(obj.time,mean(jaw{2},2),'Color',clrs.lhit,'LineWidth',2)


figure; hold on;
for j = 1:numel(params.trialid{2})
    toplot = mySmooth(jaw{1}(:,j),sm);
    patchline(obj.time,toplot,'EdgeColor',clrs.rhit,'EdgeAlpha',alph,'LineWidth',lw)
end
for j = 1:numel(params.trialid{3})
    toplot = mySmooth(jaw{2}(:,j),sm);
    patchline(obj.time,toplot,'EdgeColor',clrs.lhit,'EdgeAlpha',alph,'LineWidth',lw)
end


%% choice decoding in single trial cd projections
% cd early (use late sample for decoding)

% cd late (use all delay for decoding)

% cd go (use 0 to 0.5 sec for deocoding












