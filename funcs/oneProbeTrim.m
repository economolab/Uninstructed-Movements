function [obj,params] = oneProbeTrim(obj,params)
obj.psth = obj.psth{1};
obj.trialdat = obj.trialdat{1};
obj.presampleFR = obj.presampleFR{1};
obj.presampleSigma = obj.presampleSigma{1};
params.cluid = params.cluid{1};
end