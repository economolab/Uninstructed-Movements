function [obj,params] = loadSessionData(meta,params)

obj = loadObjs(meta);


for sessix = 1:numel(meta)
    for prbix = 1:numel(params.probe{sessix})
        disp('______________________________________________________')
        disp(['Processing data for session ' [meta(sessix).anm '_' meta(sessix).date ' | Probe' num2str(params.probe{sessix}(prbix))   ]])
        disp(' ')

        prbnum = params.probe{sessix}(prbix);

        [sessparams{sessix,prbix},sessobj{sessix,prbix}] = processData(obj(sessix),params,prbnum);
    end
end

% clean up sessparams and sessobj
for sessix = 1:numel(meta)
    params.trialid{sessix} = sessparams{sessix}.trialid;

    if numel(params.probe{sessix}) == 1
        params.cluid{sessix} = sessparams{sessix,1}.cluid{params.probe{sessix}};

        objs(sessix) = sessobj{sessix,1};
        objs(sessix).psth = objs(sessix).psth{params.probe{sessix}};
        objs(sessix).trialdat = objs(sessix).trialdat{params.probe{sessix}};
        objs(sessix).presampleFR = objs(sessix).presampleFR{params.probe{sessix}};
        objs(sessix).presampleSigma = objs(sessix).presampleSigma{params.probe{sessix}};
    elseif numel(params.probe{sessix}) == 2 % concatenate both probes worth of data

        params.cluid{sessix} = {sessparams{sessix,1}.cluid{params.probe{sessix}(1)}, sessparams{sessix,2}.cluid{params.probe{sessix}(2)} };

        objs(sessix) = sessobj{sessix,1};

        objs(sessix).psth = cat(2, objs(sessix).psth{1}, sessobj{sessix,2}.psth{2}); 
        objs(sessix).trialdat = cat(2, objs(sessix).trialdat{1}, sessobj{sessix,2}.trialdat{2}); 
        objs(sessix).presampleFR = cat(1, objs(sessix).presampleFR{1}, sessobj{sessix,2}.presampleFR{2}); 
        objs(sessix).presampleSigma = cat(1, objs(sessix).presampleSigma{1}, sessobj{sessix,2}.presampleSigma{2}); 
    end
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')


clear obj
obj = objs;
clear objs


%% convert params to a struct array (just for convenience)

temp = params;
clear params
for sessix = 1:numel(meta)
    temp2 = rmfield(temp,{'probe','trialid','cluid'});
    temp2.probe = temp.probe{sessix};
    temp2.trialid = temp.trialid{sessix};
    temp2.cluid = temp.cluid{sessix};
    params(sessix) = temp2;
end