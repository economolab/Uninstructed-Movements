function input_data = getInputData(obj,fa,data_type,np_type,varargin)

if nargin > 4
    kinfeats_reduced = varargin{1};
elseif nargin > 5
    error('too many input arguments')
end

switch np_type
    case 'optim'
        switch data_type
            case 'binned_rates'
                input_data.neural = obj.trialdat;
            case 'fa'
                input_data.neural = fa.falatents;
        end
    case 'pca'
        switch data_type
            case 'binned_rates'
                input_data.neural = obj.trialdat;
            case 'fa'
                input_data.neural = fa.falatents;
        end    
    case 'kin'
        switch data_type
            case 'binned_rates'
                input_data.neural = obj.trialdat;
                input_data.behav  = kinfeats_reduced;
            case 'fa'
                input_data.neural = fa.falatents;
                input_data.behav  = kinfeats_reduced;
        end

end
end
