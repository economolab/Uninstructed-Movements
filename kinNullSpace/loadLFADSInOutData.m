function dat = loadLFADSInOutData(pth,anm,date,runix)

pth = fullfile(pth,'lfads');

anm_date = [anm '_' date];

if isempty(runix)
    runix = getMostRecentRun(fullfile(pth,'input'),anm_date);
end

% --lfads input data
% this is a .mat file that contains the meta/params/obj used to preprocess
% lfads input data and the lfads input data itself
infn = [anm_date '_' runix '.mat'];

temp = load(fullfile(pth,'input',infn));
datin = temp.lfads; clear temp;

% --lfads output data
% there are 2 files, one for the training set, one for the validation set
% we will load both and concatenate the trials from each
outprefix = {'train',...
             'valid'};
for i = 1:numel(outprefix)
    
    outpatterns = {outprefix{i},[ anm_date '_' runix]};
    contents = dir(fullfile(pth,'output'));
    [~,mask] = patternMatchCellArray({contents.name}',outpatterns','all');
    
    if sum(mask)==1
        outfn = contents(mask).name;
    else
        error(['Found multiple matching files for lfads' outprefix{i}])
    end
    
    % the outputs of lfads are stored as hdf5 datasets
    % let's store this in a struct
    datout{i} = HDF5ToStruct(fullfile(pth,'output'),outfn);
end


% --concatenate train and valid trials after removing fields we dont need
remov = {'costs','gen_ics','gen_states','nll_bound_iwaes','nll_bound_vaes','post_g0_logvar','post_g0_mean',...
         'prior_g0_logvar','prior_g0_mean','train_steps'};
train = rmfield(datout{1},remov);
valid = rmfield(datout{2},remov);

clear datout
datout.trials = sort(cat(1,datin.train_trials,datin.valid_trials));
datout.factors = train.factors;
datout.rates = train.output_dist_params; 
% datout.factors = cat(3,train.factors,valid.factors);
% datout.rates = cat(3,train.output_dist_params,valid.output_dist_params);

% -- create one struct dat
dat = datin;
dat.trials = datout.trials;
dat.factors = permute(datout.factors,[2,1,3]);
dat.rates = permute(datout.rates,[2,1,3]);
% dat.rates = dat.rates(:,:,dat.trials); % only keep training trials (which is all trials from a session)
% dat.factors = dat.factors(:,:,dat.trials);
end