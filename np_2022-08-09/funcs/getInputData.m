function input_data = getInputData(obj,fa,gpfa,data_type,np_type,varargin)

if nargin > 5
    kinfeats_reduced = varargin{1};
elseif nargin > 6
    error('too many input arguments')
end

switch np_type
    case 'optim'
        switch data_type
            case 'binned_rates'
                input_data.neural = obj.trialdat;
            case 'fa'
                input_data.neural = fa.falatents;
            case 'gpfa'
                input_data.neural = gpfa.gpfalatents;
        end
    case 'pca'
        switch data_type
            case 'binned_rates'
                input_data.neural = obj.trialdat;
            case 'fa'
                input_data.neural = fa.falatents;
            case 'gpfa'
                input_data.neural = gpfa.gpfalatents;
        end    
    case 'kin'
        switch data_type
            case 'binned_rates'
                input_data.neural = obj.trialdat;
                input_data.behav  = kinfeats_reduced;
            case 'fa'
                input_data.neural = fa.falatents;
                input_data.behav  = kinfeats_reduced;
            case 'gpfa'
                input_data.neural = gpfa.gpfalatents;
                input_data.behav  = kinfeats_reduced;
        end

end
end
