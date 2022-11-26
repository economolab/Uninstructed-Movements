% Function for finding the average probability of jaw movement at all time
% points in a trial 

% If you input multiple sessions for an animal, it averages across all sessions 

% INPUT to 'jawkinAnimalAvg' function: (direction: 'R' or 'L', context: 'afc' for 2AFC or
% 'auw' for autowater, animal, datesToAnalyze, sessions)

% OUTPUT: saves the average probability of jaw movement at all timepoints
% in a trial.  Each column = one session


function jawdat = jawkinAnimalAvg(dir,context, animal, datesToAnalyze, sessions)
for k = 1:sessions      % For all sessions...
    
    date = datesToAnalyze(k);   % Get the current date
    pth = ['D:\JEB\DataObjects\', animal, '\Analysis\', '\data_structure_', animal,'_',date,'.mat'];
    pth = append(pth(1),pth(2),pth(3),pth(4),pth(5),pth(6),pth(7),pth(8));
    pth = string(pth);
    dat = load(pth);
    [R,L,autowater,hit,stimnum,Ntrials, traj, thisobj]  = initialize(dat);
    
    view = 2;
    posthresh = 5;      %If the y-position of the jaw crosses this threshold, the jaw is considered to have moved
    derivthresh = 0.3;  %If the velocity of the jaw crosses this threshold, the jaw is considered to be moving
    edges = 0:0.005:5.5;
    Ntottrials = Ntrials;
    medpos = zeros(Ntrials, 1);
    good_ix = NaN(Ntrials, 1);
    count = 0;
    
    %Find the trials that have numerical values for frameTimes (aka no NaNs)
    for i = 1:Ntrials
        if ~isnan(traj(i).frameTimes)
            good_ix(i,1) = i;       %If the trial is okay, note the trial number
        end
    end
    
    %Get trial numbers of usable trials
    bp = thisobj.bp;
    params.cond{1} = bp.R&bp.hit&~bp.autowater&~bp.stim.num;
    params.cond{2} = bp.L&bp.hit&~bp.autowater&~bp.stim.num;
    params.cond{3} = bp.R&bp.hit&bp.autowater&~bp.stim.num;
    params.cond{4} = bp.L&bp.hit&bp.autowater&~bp.stim.num;
    use_trials = find(~isnan(good_ix))';
    Ntrials = 0;
    ix = NaN(1,numel(good_ix));
    
    for i=1:numel(use_trials)
        if context == 'afc' % 2afc
            if dir == 'R' && params.cond{1}(use_trials(i))
                count = count+1;
                Ntrials = Ntrials +1;
                ix(count) = use_trials(i);
                clr = 'b';
            elseif dir == 'L'&& params.cond{2}(use_trials(i))
                count = count+1;
                Ntrials = Ntrials +1;
                ix(count)= use_trials(i);
                clr ='r';
            end
        else % aw
            if dir == 'R' && R(use_trials(i)) && params.cond{3}(use_trials(i))
                count = count+1;
                Ntrials = Ntrials +1;
                ix(count)= use_trials(i);
                clr = [ 0.4940    0.1840    0.5560];
                %              str = 'Right Autowater';
            elseif dir == 'L'&& L(use_trials(i)) && params.cond{4}(use_trials(i))
                count = count+1;
                Ntrials = Ntrials +1;
                ix(count)= use_trials(i);
                clr = [0.9290 0.6940 0.1250];
                %              str = 'Left Autowater';
            end
        end
    end
    ix = ix(find(~isnan(ix)));
    jaw = NaN(numel(edges), length(use_trials));
    
    for i = 1:Ntottrials
        ts = MySmooth(traj(i).ts(:, 2, 2), 21);      %Get the y-position of the jaw for the ith trial
        ts(isnan(ts)) = 0;
        medpos(i) = prctile(ts(10:end), 10);                %Finds the 10th percentile for jaw position in each trial
    end
    
    
    for i = 1:numel((ix))
        q = use_trials(find(use_trials == ix(i)));
        
        ts = MySmooth(traj(q).ts(:, 2, 2), 21);
        
        tsinterp = interp1(traj(q).frameTimes-0.5, ts, edges);   %Linear interpolation of jaw position to keep number of time points consistent across trials
        
        basederiv = nanmedian(diff(tsinterp));          %Find the median jaw velocity (aka baseline)
        
        basepos = prctile(tsinterp(10:end), 5);         %Find median jaw position
        
        %Find when the difference between the jaw velocity and the
        %baseline jaw
        %velocity is above a given threshold (when is jaw moving?)
        jaw(2:end, q) = abs(diff(tsinterp)-basederiv)>derivthresh;% | abs(tsinterp(2:end)-basepos)>posthresh;
        
    end
    
    jawdat(:,k) = nanmean(jaw,2);
    jawdat(:,k) = medfilt1(jawdat(:,k),10);
end 
y = nanmean(jawdat,2);
plot(edges, y,'color', clr, 'LineWidth', 2);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:500:2500, ...
    'LineWidth', 1)
end % jawkin
%%
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