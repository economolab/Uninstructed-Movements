function [obj,cluid] = removeLowFRClusters(obj,cluid,lowFR,prbnum)
% % Remove low-firing rate units, e.g., all those firing less than 5
%   spikes per second on average across all trials.
%   
%   The fitted observation noise (diagonal element of R) for a
%   low-firing rate unit will be small, so the neural trajectory may
%   show a deflection each time the neuron spikes.

meanFRs = nanmean(nanmean(obj.psth{prbnum},3));
use = meanFRs > lowFR;

% remove low fr clusters
cluid = cluid(use);
obj.psth{prbnum} = obj.psth{prbnum}(:,use,:);
obj.trialdat{prbnum} = obj.trialdat{prbnum}(:,use,:);

end % removeLowFRClusters