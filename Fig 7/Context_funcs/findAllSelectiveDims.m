function [selectiveDims,projdif] = findAllSelectiveDims(currez,spacename,params,modparams,conds2use)
trialproj = currez.(spacename);
    nTrials = size(trialproj,2);
    nDims = size(trialproj,3);
    PresampProj = zeros(nDims,nTrials);
    for d = 1:nDims                                               % For each dimension...
        for t = 1:nTrials                                                   % Go through all of the trials 
            PresampProj(d,t) = mean(trialproj(modparams.start:modparams.stop,t,d),1);             % Take the avg proj onto this dimension before the sample tone 
        end
    end

    temp = PresampProj;                                        % (cells x trials)
    
    % Sub-sample trials
    for c = 1:length(conds2use)
        cond = conds2use(c);                                   
        trix = params.trialid{cond};
        trix2use = randsample(trix,modparams.subTrials);
        epochAvg{c} = temp(:,trix2use);
    end
    totalprojAfc = sum(epochAvg{1},2);
    totalprojAW = sum(epochAvg{2},2);
    projdif = (abs(totalprojAfc)>abs(totalprojAW));       % 1 if AFC-preferring; 0 if AW-preferring
   
    [selectiveDims] = getSelectiveDimensions(epochAvg,modparams.sig);
