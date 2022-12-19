
function [cdContext_mode] = findCDContext(objs,sessix,ev,params)

% cdContext mode
e1 = mode(ev.sample) - 0.42 - median(ev.(params.alignEvent));
e2 = mode(ev.goCue) - 0.02 - median(ev.(params.alignEvent));

times.context = objs{sessix}.time>e1 & objs{sessix}.time<e2;
cdContext_mode = calcCD(objs{sessix},times.context);
end