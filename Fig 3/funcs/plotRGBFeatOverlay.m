function plotRGBFeatOverlay(sessix,cond2use,params,trix2use,featix,feat2use,kin,sm,ptiles,meta,times,goix)
allkin = [];
cond = cond2use;
condtrix = params(sessix).trialid{cond};
ntrials = length(condtrix);
randtrix = randsample(condtrix,trix2use);
for f = 1:length(featix)
    currfeat = featix(f);
    %         if f~=3
    currkin = mySmooth(kin(sessix).dat_std(times.startix:goix,randtrix,currfeat),sm);
    currkin = abs(currkin);
    %         else
    %             % Mean-center paw
    %             tempkin = kin(sessix).dat_std(:,randtrix,currfeat);
    %             presampME = squeeze(mean(tempkin(times.startix:times.stopix,:,:),1,'omitnan'));
    %             avgpresampME = mean(presampME,'omitnan');
    %             tempkin = tempkin-avgpresampME;
    %             currkin = mySmooth(tempkin,sm);
    %             currkin = currkin(times.startix:goix,:);
    %         end

    % Max normalize current feature
    % ^ Don't actually want to max normalize because then will be
    % normalizing by an outlier value probably
    % Want to normalize to the 90-99th percentile of values to account
    % for more of the data
    abskin = abs(currkin);
    normkin = abskin./prctile(abskin(:), ptiles(f));
    normkin(normkin>1) = 1;                                              % Will end up with values greater than 1 in this case--set these to 1
    %         normkin = 1-normkin;
    %         maxkin = max(currkin,[],"all");
    %         currkin = abs(currkin./maxkin);

    allkin = cat(3,allkin,normkin);                                      % Concatenate across features (trials x time x feat)
end

allkin = permute(allkin,[2 1 3]);                                        % (time x trials x feat/RGB)
RI = imref2d(size(allkin));
RI.XWorldLimits = [0 3];
RI.YWorldLimits = [2 5];
IMref = imshow(allkin, RI,'InitialMagnification','fit');
title(['RGB = ' feat2use '; ' meta(sessix).anm meta(sessix).date])
sgtitle(params(sessix).condition{cond2use})