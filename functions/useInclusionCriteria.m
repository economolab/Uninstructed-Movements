function [meta,objs,params] = useInclusionCriteria(objs,params,meta)
% remove sessions if:
% 1) less than 40 trials of rhit and lhit each (same as
% 2) atleast 10 clusters of quality that is not noise
% 3) Make sure traj is updated with new DLC model
use = false(size(objs));
for i = 1:numel(use)
    check(1) = numel(params.trialid{i}{1}) > 40;
    check(2) = numel(params.trialid{i}{2}) > 30;
    check(3) = numel(params.cluid{i}) >= 10;
    check(4) = length(objs{i}.traj{1}(1).featNames) == 7;
    if all(check)
        use(i) = true;
    end
end
meta = meta(use);
objs = objs(use);
params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);
end