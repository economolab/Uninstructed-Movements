function obj = findMoveOnset(obj)
% returns idxs of movement onset for an entire session

[~,movTime] = getMoveIdx(obj);

dt = 1/200;


obj.bp.ev.moveOnset = nan(size(obj.bp.ev.goCue));
for trial = 1:obj.bp.Ntrials

    movtime = movTime{trial};
    gocue = obj.bp.ev.goCue(trial);
    
    % find closest movtime to go cue time
    [~,moveonsetix] = min(abs(movtime - gocue));
    moveonset = movtime(moveonsetix);
    
    
    % then iteratively move backwards from go cue and find where animal is
    % no longer moving. move onset is this time the animal is no longer
    % moving + dt
    
    intersect = true;
    tempmoveonset = moveonset;
    while intersect % while the moveonset time and movetime vectors intersect
        tempmoveonset = tempmoveonset - dt;
        
        [~,intersectix] = min(abs(movtime - tempmoveonset));
        intersecttime = movtime(intersectix);
        timediff = abs(intersecttime - tempmoveonset);
                
        if timediff <= dt
        % if diff b/w intersecttime and tempmoveonset time is small, but tempmoveonset < 0, quit
            if tempmoveonset <= 0
                obj.bp.ev.moveOnset(trial) = gocue;
                intersect = false;
                break
            end
        % if diff b/w intersecttime and tempmoveonset time is small, keep going
            continue
        else
            % if P doesn't exist in movtime, movonset = tempmoveonset + dt
            % and intersect = false
            obj.bp.ev.moveOnset(trial) = tempmoveonset + dt;
            intersect = false;
        end
    end

end

%% if move onset time == 0, set it to go cue. These are exceptions that are hard to deal with
for i = 1:obj.bp.Ntrials
    if obj.bp.ev.moveOnset(i) < 0.01
        obj.bp.ev.moveOnset(i) = obj.bp.ev.goCue(i);
    end
end

%% validation

% dt = 1/400;
% view = 1; % side cam
% tongfeat = 1;
% jawfeat = 2; 
% sm = 7; % smoothing window
% vid = obj.traj{view};
% for j = 1:obj.bp.Ntrials
%     % get jaw trajectory
%     p = vid(j).ts(:, 3, jawfeat);
%     Nframes = size(vid(j).ts, 1);
%     dat = zeros(Nframes, 2);
%     dat(:,1) = (1:Nframes)./(1/dt); % time
%     dat(:,2) = mySmooth(vid(j).ts(:,2,jawfeat), sm); % z
%     dat(p<0.9,2) = NaN;
%     
%     % fill nans
%     for i = 2
%         dat(:,i) = fillmissing(dat(:,i),'makima');
%     end
%     
%     % mean-center and normalize jaw x and z pos
%     for i = 2
%         dat(:,i) = (dat(:,i) - nanmin(dat(:,i))) / (nanmax(dat(:,i)) - nanmin(dat(:,i)));
%         dat(:,i) = dat(:,i) - nanmean(dat(:,i));
%         dat(1:5,i) = dat(6,i);
%     end    
%     
%     % plot
%     close all
%     figure;
%     plot(dat(:,1),dat(:,2)); 
%     hold on
%     xline(obj.bp.ev.moveOnset(j),'r')
%     pause
%     
% end



end % findMoveOnset