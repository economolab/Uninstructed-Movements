function [meta,params,obj,dat] = getNeuralActivity(meta,params)

switch params.lfads_or_fa
    case 'lfads'
        % if we're using lfads output, all preprocessing has been done, so we
        % just need to load some data and move on to processing the video features
        
        % load lfads input and output data
        dat = loadLFADSInOutData(meta.datapth,meta.anm,meta.date,params.lfads_run);
        obj = dat.obj;
        % get params used to create lfads input data
        temp = dat.preprocess_params;
        temp.lfads_or_fa = params.lfads_or_fa;
        temp.lfads_run = params.lfads_run;
        params = temp; clear temp
        % get rid of stuff we no longer need
        dat = rmfield(dat,{'preprocess_params','obj','train_percent','bin_spikes',...
            'train_trials','valid_trials',...
            'train_data','valid_data'});
        
        disp('DONE LOADING LFADS DATA')
        
    case 'fa'
        % if we're performing factor analysis we will use the params and
        % session info defined in PARAMETERS section to:
        %  1) load and preprocess data
        %  2) perform factor analysis on single trials to obtain as many
        %     factors as needed to explain 75% of variance
        %  3) smooth factors and single trials
        
        % load data
        obj = loadDataForFA(meta);
        % preprocess data
        [obj,params] = processDataForFA(obj,params);
        % perform factor analysis, also trims trials we're using to those
        % in params.condition
        dat = doFA(obj,params);
        % smooth factors and single trials
        dat = smoothPostFA(dat,obj,params);
        
        disp('DONE LOADING DATA AND PERFORMING FACTOR ANALYSIS')
        
    otherwise
        error('params.lfads_or_fa should be set to either `lfads` or `fa`')
end

end