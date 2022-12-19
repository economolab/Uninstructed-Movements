% validateTimeWarp

% for each trial, plot all licks, and warped licks

view = 1; % use side cam
c = 1; % use x coord
feat = 1; % tongue in side cam

dt = 1/400;

figure;
for trix = 1:obj.bp.Ntrials
    gc = obj.bp.ev.goCue(trix);
    ls = lickStart{trix}; % current trial lick(lix) start times
    le = lickEnd{trix}; % current trial lick(lix) end times
    ld = lickDur{trix}; % current trial lick(lix) durations
    
    for lix = 1:numel(lickStart{trix}) % lick index for current trial
        
        p_cur = p{trix,lix}; % current fit parameters for current trial and lick number
        
        % median values for current lick
        mls = med.lickStart(lix); % median lick start time
        mle = med.lickEnd(lix); % median lick end time
        mld = med.lickDur(lix); % median lick duration
        
        
        % warp
        if lix==1
            x = gc;
            y = le(lix);
            warptm = polyval(p_cur,[x y]);
        else
            x = le(lix-1);
            y = le(lix);
            warptm = polyval(p_cur,[x y]);
        end

        % plot
        hold on;
        plot([x y],[0+lix 0+lix],'-k','LineWidth',3.5)
        plot([warptm(1) warptm(2)], [0+lix 0+lix], '-m','LineWidth',2)
        xline(mode(obj.bp.ev.goCue),'--g','LineWidth',2)
        xline(mle,'--r','LIneWidth',2)
    end
    xlim([2 7])
    ylim([-2 params.nLicks])
    title(['All Licks - Trial:  ' num2str(trix)])
    xlabel('Time')
    ax = gca;
    ax.FontSize = 20;
    hold off
    pause
    clf
end