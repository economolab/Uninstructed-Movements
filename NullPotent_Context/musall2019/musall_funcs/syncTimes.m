function synced_dat = syncTimes(dattm,dat,newtm)
% interp time from video time (dattm) to newtm (neural data binned time)
% dat is (time,trials,dims)

sizeVidR = size(dat);
synced_dat = zeros(numel(newtm),size(dat,2),size(dat,3));
for iDim = 1:size(dat,3)
    for iTrial = 1:size(dat,2)
        synced_dat(:,iTrial,iDim) = interp1(dattm, dat(:,iTrial,iDim), newtm);
    end
end
