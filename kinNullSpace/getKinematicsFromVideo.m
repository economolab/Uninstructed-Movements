function [kin,featLeg] = getKinematicsFromVideo(obj,params,trialnums)


feats = params.traj_features;

taxis = obj.time + params.advance_movement;


% get x/y position and velocity for each feature in feats

xpos = cell(2,1); % one entry for each view
ypos = cell(2,1); 

xvel = cell(2,1);
yvel = cell(2,1);

for viewix = 1:numel(feats) % loop through cam0 and cam1
    xpos{viewix} = nan(numel(obj.time),numel(trialnums));
    ypos{viewix} = nan(numel(obj.time),numel(trialnums));
    xvel{viewix} = nan(numel(obj.time),numel(trialnums));
    yvel{viewix} = nan(numel(obj.time),numel(trialnums));
    for featix = 1:numel(feats{viewix}) % loop through number of features for current cam
        feat = feats{viewix}{featix};
        [xpos{viewix}(:,:,featix), ypos{viewix}(:,:,featix)] = findPosition(taxis, obj, trialnums, viewix, feat, params.alignEvent);
        [xvel{viewix}(:,:,featix), yvel{viewix}(:,:,featix)] = findVelocity(taxis, obj, trialnums, viewix, feat, params.alignEvent);
    end
end

% compute displacement from origin from x/y pos for each feature

featDisp = cell(2,1);
for viewix = 1:numel(feats) % loop through cam0 and cam1
    featDisp{viewix} = nan(size(xpos{viewix}));
    for featix = 1:numel(feats{viewix}) % loop through number of features for current cam
        tempx = xpos{viewix}(:,:,featix);
        tempy = ypos{viewix}(:,:,featix);
        featDisp{viewix}(:,:,featix) = sqrt(tempx.^2 + tempy.^2);
    end
end

% concatenate features into one big matrix of size (time,features,trials)
% (displacement + xvel + yvel)

featDispMat = cat(3, featDisp{:});
xvelMat = cat(3, xvel{:});
yvelMat = cat(3, yvel{:});

kin = cat(3,featDispMat,xvelMat,yvelMat);
kin = permute(kin,[1,3,2]);

kin(1:10,:,:) = 0; % remove any smoothing artifacts at beginning of trials

% feature legend

dispFeatNames = cell(2,1);
xvelFeatNames = cell(2,1);
yvelFeatNames = cell(2,1);
for viewix = 1:numel(feats) % loop through cam0 and cam1
    for featix = 1:numel(feats{viewix}) % loop through number of features for current cam
        feat = feats{viewix}{featix};
        dispFeatNames{viewix}{featix} = [feat '_disp_' num2str(viewix-1)];
        xvelFeatNames{viewix}{featix} = [feat '_xvel_' num2str(viewix-1)];
        yvelFeatNames{viewix}{featix} = [feat '_yvel_' num2str(viewix-1)];
    end
end
disps = cat(2, dispFeatNames{:})';
xvels = cat(2, xvelFeatNames{:})';
yvels = cat(2, yvelFeatNames{:})';

featLeg = cat(1,disps,xvels,yvels);


end
