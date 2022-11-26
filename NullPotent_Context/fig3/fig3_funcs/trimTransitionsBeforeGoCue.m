function [mdat,mdat_leg,qdat,qdat_leg] = trimTransitionsBeforeGoCue(m2q,q2m,obj,params)

for sessix = 1:numel(m2q)
    align = obj(sessix).bp.ev.(params(sessix).alignEvent);
    for trix = 1:numel(m2q(sessix).moveStart)
        if isempty(m2q(sessix).moveStart{trix})
            continue;
        end

        [~,gcix] = min(abs(obj(sessix).time)); % time aligned to go cue
        
        % each col is a m2q bout within a trial, rows are moveStart, moveEnd, etc.
        dat(1,:) = m2q(sessix).moveStart{trix};
        dat(2,:) = m2q(sessix).moveEnd{trix};
        dat(3,:) = m2q(sessix).quietStart{trix};
        dat(4,:) = m2q(sessix).quietEnd{trix};
        
        % only use m2q bouts before the go cue
        use = all(dat < gcix);
        
        mdat{sessix}{trix} = dat(:,use);
        clear dat

    end
end

for sessix = 1:numel(q2m)
    align = obj(sessix).bp.ev.(params(sessix).alignEvent);
    for trix = 1:numel(q2m(sessix).moveStart)
        if isempty(q2m(sessix).quietStart{trix})
            continue;
        end

        [~,gcix] = min(abs(obj(sessix).time)); % time aligned to go cue
        
        % each col is a q2m bout within a trial, rows are moveStart, moveEnd, etc.
        dat(1,:) = q2m(sessix).quietStart{trix};
        dat(2,:) = q2m(sessix).quietEnd{trix};
        dat(3,:) = q2m(sessix).moveStart{trix};
        dat(4,:) = q2m(sessix).moveEnd{trix};
        
        % only use q2m bouts before the go cue
        use = all(dat < gcix);
        
        qdat{sessix}{trix} = dat(:,use);
        clear dat

        % trim quiet-to-move transitions to size (each epoch tbout long)
%         qdat = q2m(sessix);
    end
end

mdat_leg = {'movestart';'moveend';'quietstart';'quietend'};
qdat_leg = {'quietstart';'quietend';'movestart';'moveend'};



end