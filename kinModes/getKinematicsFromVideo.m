function kin = getKinematicsFromVideo(obj,params,trialnums)


feats = params.traj_features;

taxis = obj.time;

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

% TODO
% handle nans (it will be different for tongue and other features)
% concatenate all features into a big matrix
% - should we be calculating displacement from origin and get rid of
%    xpos/ypos (would just have displacement for each feature)
% - should we be calculating speed from x/yvel and just have one var per
%    feature
% make a legend for the feature entries of the big matrix (time,features,trials)

kin = [];
featLeg = {};

end
