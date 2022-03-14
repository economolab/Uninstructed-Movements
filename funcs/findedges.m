function edges = findedges(time,bp,dt,epoch,trial,alignEvent)
    % find histogram bin edges for a specific trial and epoch
    % used to find activity modes / coding vectors when trial and/or epoch
    % lengths differ between trials
    switch epoch
        case 'presample'
            e1 = bp.ev.sample(trial) - 0.5; 
            e2 = bp.ev.sample(trial) - 0.05;
        case 'sample'
            e1 = bp.ev.sample(trial);
            e2 = bp.ev.delay(trial) - 0.3;
        case 'delay'
            e1 = bp.ev.delay(trial);
            e2 = bp.ev.goCue(trial) - 0.05;
        case 'prego'
            e1 = bp.ev.goCue(trial) - 0.1;
            e2 = bp.ev.goCue(trial); 
        case 'postgo'
            e1 = bp.ev.goCue(trial);
            e2 = bp.ev.goCue(trial) + 0.1; 
        case 'outcome'
            e1 = bp.ev.goCue(trial);
            e2 = bp.ev.goCue(trial) + 1.3;
        case 'moveonset'
            e1 = bp.ev.moveOnset(trial);
            e2 = bp.ev.moveOnset(trial) + 0.2;
        case 'action'
            e1 = bp.ev.goCue(trial) + 0.1;
            e2 = bp.ev.goCue(trial) + 0.3;
    end
   % align edges to params.alignEvent
   e1 = e1 - bp.ev.(alignEvent)(trial);
   e2 = e2 - bp.ev.(alignEvent)(trial);
   % define edges based on e1 and e2 for the current trial
   [~,e1ix] = min(abs(time - e1));
   [~,e2ix] = min(abs(time - e2));
   edges = [e1ix e2ix]; 
%    edges = e1:dt:e2;
   
end % findedges