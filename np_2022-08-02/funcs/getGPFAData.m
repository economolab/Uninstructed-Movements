function gpfa = getGPFAData(meta,runix)

pth = meta.datapth;
anm = meta.anm;
date = meta.date;

pth = fullfile(pth,'gpfa');

anm_date = [anm '_' date];

if isempty(runix)
    runix = getMostRecentRun(fullfile(pth,'input'),anm_date);
end

infn = [anm_date '_' runix '.mat'];

gpfa = load(fullfile(pth,infn));



end