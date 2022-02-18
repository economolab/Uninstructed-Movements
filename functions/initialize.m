function [R,L,autowater,hit,stimnum,Ntrials, traj, thisobj]  = initialize(dat)

thisobj = dat.obj;
R = thisobj.bp.R;
L = thisobj.bp.L;
autowater = thisobj.bp.autowater;
hit = thisobj.bp.hit;
early = thisobj.bp.early;
lickL = thisobj.bp.ev.lickL;
lickR = thisobj.bp.ev.lickR;
stimnum = thisobj.bp.stim.num;
no = thisobj.bp.no;
traj = thisobj.traj{1};

Ntrials = thisobj.bp.Ntrials;
end % initialize