clear,clc,close all
addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMALs/SESSIONs TO LOAD
meta(1).anm = 'JEB7';
meta(1).date = '2021-04-29';

meta(end+1).anm = 'JEB7';
meta(end).date = '2021-04-30';

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

% % meta(end+1).anm = 'EKH3';
% % meta(end).date = '2021-08-07';

meta(end+1).anm = 'EKH1';
meta(end).date = '2021-08-07';

meta = assignDataPath(meta);

% ------------------------------------------------------------------ %
% get some parameters that are common to all sessions/analyses
dfparams = getDefaultParams();

%% NEURAL ACTIVITY

% data objects contain
% 1) processed obj
% 2) processing params
% 3) factor analysis results
% 4) params used for factor analysis
% 5) meta data that we don't really need

[obj,fa,params,me,kin,kinfeats,kinfeats_reduced] = loadProcessedData(meta);

%% analysis type
np_type = 'optim_psth_me'; % 'optim'  'kin'  'pca'  'optim_psth'   'optim_psth_me'
data_type = 'binned_rates';  % 'fa'   'binned_rates'

% optim: elsayed method on single trials with move and non-move times
% kin:   kaufman method on single trials with move and non-move times
% pca:   nullPCs estimated from non-move times, potentPCs found from
%        residual neural activity after removing nullPCs



%% KINEMATICS

% if strcmpi(np_type,'kin') && exist('kinfeats','var') == 0
%     [kin,kinfeats,kinfeats_reduced] = getKin(meta,obj,dfparams,params,me);
% end


%% null/potent

clearvars -except meta params obj fa dfparams me kin kinfeats kinfeats_reduced np_type data_type

switch np_type

    case 'optim'
        for i = 1:numel(meta)
            input_data = getInputData(obj(i),fa(i),data_type,np_type);
            input_data = input_data.neural;
            rez(i) = elsayedNullandPotentSpace(obj(i),input_data,me(i),params(i));
        end

    case 'kin'
        for i = 1:numel(meta)
            input_data = getInputData(obj(i),fa(i),data_type,np_type,kinfeats_reduced{i});
            input_data.factors = input_data.neural;
            input_data.feats   = input_data.behav;
            rez(i) = kinematicNullAndPotentSpace(obj(i),input_data,dfparams,params(i),me(i));
        end

    case 'pca'
        for i = 1:numel(meta)
            input_data = getInputData(obj(i),fa(i),data_type,np_type);
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



% %% variance explained plots
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

    plotCD_psth(rez,obj,params,times);

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
        plotSessionCD(rez(i),obj(i),params(i),times)
%         % null space cds
%         plotSessionNullSpaceCD(rez(i),obj(i),params(i),times)
%         % potent space cds
%         plotSessionPotentSpaceCD(rez(i),obj(i),params(i),times)
    end

end

%% TODO

% 1. plot cds for each session separately as well
% - mean center null/potent_psths before calculating CDs
% 2. missing ve code in activity modes calc
% 3. kaufman method (full rank??), no null space
% 4. latency b/w neural activity and motion energy
% - use binned firing rates instead
% no selectivity in response period in n/p projections?
% selectivity before sample period in n/p projections??
% reproduce psth version with our data
% cd's for optim_psth - now that I understand this method better, it kinda
% seems like it's a good method, since on trial averaged data, we estimate
%       what movement looks like in the response period and then project data
%       onto that. 
% cluster trials based on motion energy, use those for optim_psth_me
% when combingin sessions for CDs, make sure one trial type is higher
% dont use normalizeInRange when plotting CDs





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
sessix = 2;
dim = 1;
plotME_PotentDim(obj,me,rez,sessix,dim,data_type)
set(gca,'XLim',[256.3261  264.4115])


% % same for a null dimension
% sessix = 2;
% dim = 1;
% plotME_NullDim(obj,me,rez,sessix,dim,data_type)
% set(gca,'XLim',[256.3261  264.4115])


















