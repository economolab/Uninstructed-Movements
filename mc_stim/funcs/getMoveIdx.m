function [movix,movtime] = getMoveIdx(obj,params)
% returns cell array mov where each entry is an (Nframesx1) logical vector
% labeling where in a trial the animal is moving

% algorithm:
% 1) set timepoint to true if tongue is visible
% 2) set timepoint to true if jaw speed exceeds a threshold
% 3) set timepoint true if variance in jaw position exceeds a threshold

% movIdx is 1) or 2) or 3) 

%% params
% 
% figure(1)
% set(gcf,'units','normalized','Position',[0.2 0.1 0.7 0.8]); 


dt = 1/400;
view = 1; % side cam
tongfeat = 1;
jawfeat = 2; 
sm = 7; % smoothing window

jawSpeedPercent = 0.35; % defines jaw speed threshold above which is considered movement
jawPosVarPercent = 0.05; % defines variance in jaw position thresh, above which is considered movement



vid = obj.traj{view};
movIdx = cell(obj.bp.Ntrials,1);
movTime = movIdx;
for j = 1:obj.bp.Ntrials
    %% get tongue trajectory
    
    p = vid(j).ts(:, 3, tongfeat);
    
    Nframes = size(vid(j).ts, 1);
    dat = zeros(Nframes, 2);
    dat(:,1) = (1:Nframes)./(1/dt); % time
    dat(:,2) = mySmooth(vid(j).ts(:,1,tongfeat), sm); % x
    dat(:,3) = mySmooth(vid(j).ts(:,2,tongfeat), sm); % z
    dat(p<0.9,[2,3]) = NaN;
    
    %% find where tongue is visible 
    
    tongmov = false(Nframes,1);
    tongmov(~isnan(dat(:,2))) = true;
    tongmov(~isnan(dat(:,3))) = true;
    
    %% get jaw trajectory

    p = vid(j).ts(:, 3, jawfeat);
    
    Nframes = size(vid(j).ts, 1);
    dat = zeros(Nframes, 2);
    dat(:,1) = (1:Nframes) .* dt; % time
    dat(:,2) = mySmooth(vid(j).ts(:,1,jawfeat), sm); % x
    dat(:,3) = mySmooth(vid(j).ts(:,2,jawfeat), sm); % z
    dat(p<0.9,[2,3]) = NaN;
    
    %% fill nans
    
    for i = 2:3
        dat(:,i) = fillmissing(dat(:,i),'makima');
    end
    
    %% mean-center and normalize jaw x and z pos
    
    for i = 2:3
        dat(:,i) = (dat(:,i) - nanmin(dat(:,i))) / (nanmax(dat(:,i)) - nanmin(dat(:,i)));
        dat(:,i) = dat(:,i) - nanmean(dat(:,i));
        dat(1:5,i) = dat(6,i);
    end    
    
    
    %% get jaw speed
    
    vx = gradient(dat(:,2), dt);
    vy = gradient(dat(:,3), dt);
    dat(:,4) = sqrt(vx.^2 + vy.^2);
%     dat(:,4) = gradient(dat(:,3), dt);

    dat(1:50,4) = dat(51,4); 
    dat(:,4) = mySmooth(abs(dat(:,4)),200);
    
    dat(:,4) = (dat(:,4) - nanmin(dat(:,4))) / (nanmax(dat(:,4)) - nanmin(dat(:,4)));
%     dat(:,4) = dat(:,4) - nanmean(dat(:,4));
    dat(1:50,4) = dat(51,4); 

    %% find where jaw speed exceeds a threshold
    
    per1 = jawSpeedPercent;
    if min(dat(:,4)) < 0
        jawspeedthresh = per1*(max(dat(:,4)) + min(dat(:,4)));
    else
        jawspeedthresh = per1*(max(dat(:,4)) - min(dat(:,4)));
    end
    jawspeedmov = dat(:,4) > jawspeedthresh;
        
    
    %% find where variance in jaw z pos exceeds a threshold
    
    val = zeros(Nframes, 1);
    binsize = 0.5; % in seconds
    for ii = (1/dt * binsize):Nframes % for each half second bin
        val(ii) = nanvar(dat(ii-(1/dt*binsize-1):ii,3)); % variance in jaw z pos in bin
    end
    
    per2 = jawPosVarPercent;
    jawzthresh = per2*max(val);
    jawzmov = val > jawzthresh;
    
    %% get mov
    
    movIdx{j} = tongmov | jawspeedmov | jawzmov;
    movTime{j} = dat(movIdx{j},1);


     %% align to same event as neural data and only keep what's in obj.time
    
    movtime{j} = movTime{j} - obj.bp.ev.(params.alignEvent)(j);
    movtime{j} = movtime{j}((movtime{j}>=min(obj.time) & movtime{j}<=max(obj.time)));
    
    %% movIdx is then the idx in obj.time where value of obj.time is closest to movTime
    movix{j} = nan(size(movtime{j}));
    for i = 1:numel(movtime{j})
        mtime = movtime{j}(i);
        [~,movix{j}(i)] = min(abs(obj.time - mtime));
    end
    movix{j} = unique(movix{j});
    
%     %% plot
% % %     close all
%     subplot(3,1,1)
%     plot(dat(:,1),dat(:,4),'LineWidth',2.5); hold on;
%     plot(dat(jawspeedmov,1),dat(jawspeedmov,4),'.','MarkerSize',1); hold off
%     yline(jawspeedthresh,'r');
%     ylabel('jaw speed')
%     ax = gca;
%     ax.FontSize = 20;
%     
%     subplot(3,1,2)
%     title('variance in jaw z pos')
%     plot(dat(:,1),val,'LineWidth',2.5); hold on
%     plot(dat(jawzmov,1),val(jawzmov),'.','MarkerSize',1); hold off
%     yline(jawzthresh,'r');
%     ylabel('variance in jaw z pos')
%     ax = gca;
%     ax.FontSize = 20;
%     
%     subplot(3,1,3)
%     plot(dat(:,1),dat(:,3),'LineWidth',2.5); hold on
%     plot(dat(movIdx{j},1),dat(movIdx{j},3),'.','MarkerSize',1); hold off;
%     ylabel('jaw z position')
%     xlabel('time')
%     sgtitle(['Trial ' num2str(j)])
%     ax = gca;
%     ax.FontSize = 20;
%     
%     pause    
     
    
end

end % getMoveidx




















