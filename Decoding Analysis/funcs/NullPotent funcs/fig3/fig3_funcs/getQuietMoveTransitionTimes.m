function [q2m, m2q] = getQuietMoveTransitionTimes(me,params,obj,tbout)

for sessix = 1:numel(me)
    dt = params(sessix).dt;
    ntime = numel(obj(sessix).time);
    for trix = 1:size(me(sessix).data,2)
        ct1 = 1; % keep track of each valid move bout with valid pre quiet bout
        ct2 = 1; % keep track of each valid move bout with valid post quiet bout


        % find bouts of movement
        [mstart, mend, mlen, mcounts] = ZeroOnesCount(me(sessix).move(:,trix));
        
        % find all bouts of movement equal to or longer than tbout
        tlen = mlen .* dt;

        % for each movement bout, see if its longer than tbout, 
        % then if there is a non-move bout
        % either before or after it that is equal to or longer than tbout
        for i = 1:mcounts

            % check if current move bout is longer than tbout
            if ~(tlen(i) >= tbout)
                continue
            end

            % check if there's a quiet tbout before move bout
            % basically, how long to get to previous move bout, unless i==1
            if i == 1
                % how long to first move bout
                cur = mstart(i) .* dt;
                qlen = cur;
            else
                % how long to previous move bout
                cur = mstart(i) .* dt;
                prev = mend(i-1) .* dt;
                qlen = cur - prev;
            end
            if qlen >= tbout
                q2m(sessix).quietStart{trix}(ct1) = mstart(i) - (tbout/dt);
                q2m(sessix).quietEnd{trix}(ct1)   = mend(i) - 1;
                q2m(sessix).moveStart{trix}(ct1)  = mstart(i);
                q2m(sessix).moveEnd{trix}(ct1)    = mend(i);
                ct1 = ct1 + 1;
            end

            % check if there's a quiet tbout after move bout
            % basically, how long to get to next move bout, unless i==end
            if i == mcounts
                % how long to end of time
                cur = mend(i);
                qlen = (ntime - cur) .* dt;
            else
                % how long to next move bout
                cur = mend(i) .* dt;
                next = mstart(i+1) .* dt;
                qlen = next - cur;
            end
            if qlen >= tbout
                m2q(sessix).moveStart{trix}(ct2) = mstart(i);
                m2q(sessix).moveEnd{trix}(ct2)    = mend(i);
                m2q(sessix).quietStart{trix}(ct2) = mend(i) + 1;
                if i == mcounts
                    m2q(sessix).quietEnd{trix}(ct2)   = mend(i) + (tbout/dt);
                else
                    m2q(sessix).quietEnd{trix}(ct2)   = mstart(i+1) - 1;
                end
                ct2 = ct2 + 1;
            end
        end

    end
end



end