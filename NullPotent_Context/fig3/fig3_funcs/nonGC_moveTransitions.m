function [mdat,mdat_leg,qdat,qdat_leg,me] = nonGC_moveTransitions(obj,me,params,stitch_dist,purge_dist,tbout)

% 'denoise' movement bouts and purge any bouts shorter than purge_dist

newme = stitchAndPurgeMoveBouts(me,params,stitch_dist,purge_dist);

% from 'denoised' move bouts, find movement and non-movement transitions
[q2m, m2q] = getQuietMoveTransitionTimes(newme,params,obj,tbout);
% q2m/m2q : struct array (one entry per session). Each field has a cell per
% trial, and each trial is an array of time indices for move and quiet
% bout start and end positions
% q2m - quiet to move transitions
% m2q - move to quiet transitions

% trim to just pre-gocue transitions
[mdat,mdat_leg,qdat,qdat_leg] = trimTransitionsBeforeGoCue(m2q,q2m,obj,params);

me = newme; % me.move is different after stitch/purge, all else same



end % nonGC_moveTransitions




















