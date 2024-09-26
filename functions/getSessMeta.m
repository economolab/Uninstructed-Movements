function [anm,date,probenum,taxis] = getSessMeta(obj,met)
anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used
taxis = obj.time;
end