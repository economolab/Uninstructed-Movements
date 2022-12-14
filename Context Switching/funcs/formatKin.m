function reformattedKin = formatKin(kinData,params,time,cond2use,feats2use)
% Function for formatting kinematic data into same format as condition-averaged PSTHs
% End up with 'reformattedKin' which will be the average kinematic data over trials in each condition for each feature 
reformattedKin = NaN(length(time),length(feats2use),length(cond2use));     % Want kinematic data to be in format (time x nMoveFeat x conditions)

for c = 1:length(cond2use)                                                 % For each condition...
    condtrix = params.trialid{cond2use(c)};                      % Get the trials which correspond to that condition
    condkin = kinData(:,condtrix,feats2use);                     % Grab the kinematic data for desired features from those trials
    condkin = squeeze(mean(condkin,2,'omitnan'));                % Average the kinematic data across all trials from this condition
    reformattedKin(:,:,c) = condkin;                             % Store kin data in proper format
end