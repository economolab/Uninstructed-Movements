function testsplit = getTestTrials(params,cond2use,nSplits)
for sessix = 1:length(params)
% Get number of test trials that should be in each split
    for c = 1:length(cond2use)
        if c==1
            cond = 'afc';
        else
            cond = 'aw';
        end
        allcondtrix = params(sessix).trialid{cond2use(c)};
        trixPersplit = floor(length(allcondtrix)/nSplits);
        usedtrix = false(1,length(allcondtrix));
        for ss = 1:nSplits
            usabletrix = allcondtrix(~usedtrix);
            test = randsample(usabletrix,trixPersplit);
            testix = find(ismember(allcondtrix,test));
            usedtrix(testix) = 1;
            
            fn = ['Split' num2str(ss)];
            testsplit(sessix).(fn).(cond) = test;
        end
    end
end