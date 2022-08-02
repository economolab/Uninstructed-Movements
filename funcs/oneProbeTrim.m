function [obj,params] = oneProbeTrim(obj,params)
obj.psth = obj.psth{params.probe};
obj.trialdat = obj.trialdat{params.probe};
obj.presampleFR = obj.presampleFR{params.probe};
obj.presampleSigma = obj.presampleSigma{params.probe};
params.cluid = params.cluid{params.probe};

if isfield(obj,'trialspikes')
    obj.trialspikes = obj.trialspikes{params.probe};
end

end