function [kin,featLeg] = getKinematicsFromVideo(obj,params)

% ----------------------------------------------
% -- parameters --
% ----------------------------------------------
NvarsPerFeat = 4; % (xdisp,ydisp,xvel,yvel)

feats = params.traj_features; % features to extract kinematics for

taxis = obj.time + params.advance_movement; % time axis to align to

nTrials = obj.bp.Ntrials; 

% ----------------------------------------------
% -- get kinematics --
% ----------------------------------------------

xpos = cell(2,1); % one entry for each view
ypos = cell(2,1); 

xvel = cell(2,1);
yvel = cell(2,1);

for viewix = 1:numel(feats) % loop through cam0 and cam1
    xpos{viewix} = nan(numel(obj.time),nTrials);
    ypos{viewix} = nan(numel(obj.time),nTrials);
    xvel{viewix} = nan(numel(obj.time),nTrials);
    yvel{viewix} = nan(numel(obj.time),nTrials);
    for featix = 1:numel(feats{viewix}) % loop through number of features for current cam
        feat = feats{viewix}{featix};

        [xpos{viewix}(:,:,featix), ypos{viewix}(:,:,featix)] = findPosition(taxis, obj, nTrials, viewix, feat, params.alignEvent);
        [xvel{viewix}(:,:,featix), yvel{viewix}(:,:,featix)] = findVelocity(xpos{viewix}(:,:,featix), ypos{viewix}(:,:,featix), feat);
    end
end

% compute displacement from origin from x/y pos for each feature

featMat = cell(2,1);
viewNum = cell(2,1);
for viewix = 1:numel(feats) % loop through cam0 and cam1
    featMat{viewix} = nan(size(xpos{viewix}, 1), size(xpos{viewix}, 2), size(xpos{viewix}, 3).*NvarsPerFeat);
    for featix = 1:numel(feats{viewix}) % loop through number of features for current cam
        tempx = xpos{viewix}(:,:,featix);
        tempy = ypos{viewix}(:,:,featix);

%         featDisp{viewix}(:,:,featix) = sqrt((tempx.^2 + tempy.^2);
        featMat{viewix}(:,:,(featix*NvarsPerFeat)-3) = tempx;
        featMat{viewix}(:,:,featix*NvarsPerFeat-2) = tempy;
        featMat{viewix}(:,:,featix*NvarsPerFeat-1) = xvel{viewix}(:, :, featix);
        featMat{viewix}(:,:,featix*NvarsPerFeat) = yvel{viewix}(:, :, featix);
        viewNum{viewix}(featix) = viewix;
    end
end

% concatenate features into one big matrix of size (time,features*2,trials)
% (displacement + xvel + yvel)

kin = cat(3, featMat{:});
v = cat(2, viewNum{:});
kin(1:5, :, :) = kin(1:5, :, :).*0; % remove any smoothing artifacts at beginning of trials

allfeats = cat(2, feats{:});
featLeg = cell(size(kin, 3), 1);
for i = 1:numel(allfeats)
    featLeg{i*NvarsPerFeat-3} = [allfeats{i} '_xdisp_view' num2str(v(i))];
    featLeg{i*NvarsPerFeat-2} = [allfeats{i} '_ydisp_view' num2str(v(i))];
    featLeg{i*NvarsPerFeat-1} = [allfeats{i} '_xvel_view' num2str(v(i))];
    featLeg{i*NvarsPerFeat} = [allfeats{i} '_yvel_view' num2str(v(i))];
end

end
