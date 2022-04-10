function jawAngle = getJawAngle(taxis, obj, trialnums, alignEv)

%jaw angle

Ntrials = numel(trialnums);
jawAngle = zeros(numel(taxis), Ntrials);

view = 2; % bottom cam

jawIx = findDLCFeatIndex(obj.traj,view,{'jaw'});
topTongIx = findDLCFeatIndex(obj.traj,view,{'top_tongue'});
botTongIx = findDLCFeatIndex(obj.traj,view,{'bottom_tongue'});

for i = 1:Ntrials
    trix = trialnums(i);
    dat = medfilt1(obj.traj{view}(trix).ts, 3, [], 1);
    
    % find index corresponding to jaw for dat
    
    %Find the median x and y jaw position for the trial i
    jx = nanmedian(dat(:, 1, jawIx));
    jy = nanmedian(dat(:, 2, jawIx));
    
    %Find the x and y tongue tip position for the all time points in trial i
    tx = (dat(:,1,topTongIx)+dat(:,1,botTongIx))./2;    %Average between x position of top tongue and bottom tongue
    ty = (dat(:,2,topTongIx)+dat(:,2,botTongIx))./2;    %Average between y position of top tongue and bottom tongue
    
    dx = tx-jx;                             %Distance in x coordinates between tongue tip and jaw
    dy = ty-jy;                             %Distance in y coordinates between tongue tip and jaw
    len{i} = sqrt(dx.^2 + dy.^2);           %Length of tongue for all points in trial i
    ang{i} = atan(dy./dx);                  %Angle of tongue for all points in trial i
    
    ang{i}(dx<0 & dy>0) = ang{i}(dx<0 & dy>0) + pi;     %Correction for the quadrant that the angle lies in
    ang{i}(dx<0 & dy<0) = ang{i}(dx<0 & dy<0) - pi;
    
    jawAngle(:, i) = interp1(obj.traj{view}(trix).frameTimes-0.5-obj.bp.ev.(alignEv)(trix), ang{i}, taxis);
end







