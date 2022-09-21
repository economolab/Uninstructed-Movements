% Munib's latest version for finding CDearly, late, and go
function [cdEarly_mode, cdLate_mode, cdGo_mode] = findAllCDs(rez,sessix,ev,params)

% cd early mode
e1 = mode(ev.sample) + 0.4 - mode(ev.(params.alignEvent));
e2 = mode(ev.sample) + 0.8 - mode(ev.(params.alignEvent));

times.early = rez(sessix).time>e1 & rez(sessix).time<e2;
cdEarly_mode = calcCD(rez(sessix),times.early);


% cd late mode
e1 = mode(ev.goCue) - 0.5 - mode(ev.(params.alignEvent));
e2 = mode(ev.goCue) - 0.1 - mode(ev.(params.alignEvent));

times.late = rez(sessix).time>e1 & rez(sessix).time<e2;
cdLate_mode = calcCD(rez(sessix),times.late);


% cd go mode
e1 = mode(ev.goCue) + 0.02 - mode(ev.(params.alignEvent));
e2 = mode(ev.goCue) + 0.42 - mode(ev.(params.alignEvent));

times.go = rez(sessix).time>e1 & rez(sessix).time<e2;
cdGo_mode = calcCD(rez(sessix),times.go);
end