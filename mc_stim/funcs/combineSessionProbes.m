% INPUTS: meta, objs, anm (the animal name [as a string] that you want to
% want to combine probe data for), dates (the dates [as a cell of strings]
% that you want to combine probe data for)

function [meta, objs,params] = combineSessionProbes(meta,objs,params,anm,dates)
alldates = cell(1,length(meta));  
allanm = cell(1,length(meta));
for m = 1:length(meta)                          % Get all of the dates and animal names listed in meta
    met = meta(m);
    alldates{m} = met.date;
    allanm{m} = met.anm;
end

% Concatenate cells/PSTHs from both probes from a session
touse = true(size(meta));
for d = 1:numel(dates)                          % For each of the dates that you specified...
    currdate = dates{d};
    ix = find(strcmp(alldates,currdate)&strcmp(allanm,anm));       % Find the meta entries that correspond to that date     
    met1 = meta(ix(1)); met2 = meta(ix(2));                        % Get the meta and objs for each probe from that date
    obj1 = objs{ix(1)}; obj2 = objs{ix(2)};
    
    obj1.psth = cat(2,obj1.psth,obj2.psth);                        % Concatenate all of the probe information for the given session
    obj1.trialdat = cat(2,obj1.trialdat,obj2.trialdat);
    obj1.presampleFR = cat(1,obj1.presampleFR,obj2.presampleFR);
    obj1.presampleSigma = cat(1,obj1.presampleSigma,obj2.presampleSigma);

    objs{ix(1)} = obj1;                         % Store this concatenated info in the first obj entry for the session
    touse(ix(1)) = 1; touse(ix(2)) = 0;         % Only keep meta and obj entries for the first probe for the session
end
meta = meta(touse);
objs = objs(touse);
params.trialid = params.trialid(touse);
params.cluid = params.cluid(touse);
end