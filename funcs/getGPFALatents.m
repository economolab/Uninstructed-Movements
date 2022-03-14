function latents = getGPFALatents(obj,meta,params)

dat = formatDat(obj);
dat = preprocessDat(dat,obj,meta.trialNum);

% GPFA
resultsPth = 'C:\Code\subspace-id\gpfa\mat_results\';
runIdx = getNextRunIdx(resultsPth); 
method = 'gpfa';
xDim = 8; 
kernSD = 30; % 30 (default) - for two-stage methods
segLength = 20; % 20 (default) - didn't notice differences in changing this
binWidth = 20; % 20 (default) - differences in changing this param depends on dt
parallelize = false; % have to fix saving stuff in neuralTraj if true
numFolds = 0;

result = neuralTraj(runIdx, dat, 'datFormat', 'seq',...
                    'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD, 'segLength', segLength,...
                    'binWidth', binWidth,'parallelize',parallelize, ...
                    'numFolds', numFolds);

[~, seqTrain] = postprocess(result);

latents = nan(seqTrain(1).T,numel(seqTrain),xDim); % (time,trials,dims)
for n = 1:numel(seqTrain) % for every trial
    for k = 1:xDim % for every latent dim
        latents(:,n,k) = seqTrain(n).xorth(k,:);
    end
end

% MAKE SOME BETTER PLOTTING FUNCTIONS (A GUI LIKE THE OTHERS WOULD BE COOL)

redTrials = meta.trialNum{1};

nPlotMax = 100;
plotEachDimVsTime(seqTrain, 'xorth', obj.time,'nPlotMax',nPlotMax,...
                  'redTrials',redTrials);
nPlotMax3D = 20;
plot3D(seqTrain, 'xorth', 'dimsToPlot', [1 2 4], 'nPlotMax',nPlotMax3D,...
                  'redTrials',redTrials);


end % getGPFALatents


function dat = formatDat(obj)
% format data for gpfa
% nth entry has fields
%                       trialId      -- unique trial identifier  
%                       T (1 x 1)    -- number of timesteps
%                       y (yDim x T) -- continuous valued data 
%                                       (Eg: binned spike counts)
for i = 1:size(obj.trialpsth,3)
    dat(i).trialId = obj.trialId(i);
    dat(i).T       = numel(obj.time);
    dat(i).y = obj.trialpsth(:,:,i)';
end
end % formatDat

function runIdx = getNextRunIdx(resultsPth)
    contents = dir(resultsPth);
    contents = contents([contents.isdir]');
    dirs = {contents.name}';
    
    runIdx = 0;
    for i = 1:numel(dirs)
        if ~contains(dirs{i},'run')
            continue
        end
        curRunIdx = str2double(dirs{i}(4:end));
        if curRunIdx > runIdx
            runIdx = curRunIdx;
        end
    end
    runIdx = runIdx + 1;
end % getNextRunIdx

function dat = preprocessDat(dat,obj,trials)
% Square-root-transform spike counts to stabilize 
% the variances (Yu et al., 2009). 

for i = 1:numel(dat)
    dat(i).y = sqrt(dat(i).y);
end

% Subtract the average across trials for each condition and time bin. 
% This procedure removes systematic contributions by the stimulus, and thus, 
% the model explains only the trial-to-trial variability
means = nan(size(obj.trialpsth,1),size(obj.trialpsth,2),numel(trials));
for i = 1:numel(trials) % each condition
    % find trial indices
    trix = ismember(obj.trialId,trials{i});
    % get means across these trials
    means(:,:,i) = mean(obj.trialpsth(:,:,trix),3);
end

for i = 1:numel(dat)
    for j = 1:numel(trials)
        if ismember(dat(i).trialId,trials{j})
            cond = j;
            break
        end
    end
    dat(i).y = dat(i).y - squeeze(means(:,:,cond))';
end

% Note that in this case both the network state x and the 
% observed (transformed) spike counts y have zero mean (over trials) in each bin
% CHECK:
% ntrials = numel(dat);
% ix = 300;
% clu = 1;
% rsum = 0;
% for i = 1:ntrials
%     rsum = rsum + dat(i).y(clu,ix); 
% end
% checkmean = rsum / ntrials


end % preprocessDat







