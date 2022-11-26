function jawAngle = getJawAngle(taxis, obj, met)

%jaw angle
trialnums= [met.trialid{1}; met.trialid{2}];

Ntrials = numel(trialnums);
jawAngle = zeros(numel(taxis), Ntrials);

for i = 1:Ntrials
    num = trialnums(i);
    dat = medfilt1(obj.traj{2}(num).ts, 3, [], 1);
    
    %Find the median x and y jaw position for the trial i
    jx = nanmedian(dat(:, 1, 8));
    jy = nanmedian(dat(:, 2, 8));
    
    %Find the x and y tongue tip position for the all time points in trial i
    tx = (dat(:, 1, 1)+dat(:, 1, 3))./2;    %Average between x position of top tongue and bottom tongue
    ty = (dat(:, 2, 1)+dat(:, 2, 3))./2;    %Average between y position of top tongue and bottom tongue
    
    dx = tx-jx;                             %Distance in x coordinates between tongue tip and jaw
    dy = ty-jy;                             %Distance in y coordinates between tongue tip and jaw
    len{i} = sqrt(dx.^2 + dy.^2);           %Length of tongue for all points in trial i
    ang{i} = atan(dy./dx);                  %Angle of tongue for all points in trial i
    
    ang{i}(dx<0 & dy>0) = ang{i}(dx<0 & dy>0) + pi;     %Correction for the quadrant that the angle lies in
    ang{i}(dx<0 & dy<0) = ang{i}(dx<0 & dy<0) - pi;
    
    jawAngle(:, i) = interp1(obj.traj{2}(num).frameTimes-0.5-mode(obj.bp.ev.goCue), ang{i}, taxis);
end