function [kin,kinfeats] = getKin(meta,obj,dfparams,params)

% getKinematicsFromVideo() returns 2 variables
% - kin: feature matrix of size (time,trials,features). Features defined in
%         params.traj_features
% - featLeg: legend corresponding to features in kin (for 2nd dimension)
for i = 1:numel(meta)
    [kin(i).dat,kin(i).featLeg] = getKinematicsFromVideo(obj(i),dfparams,params(i),1:obj(i).bp.Ntrials);


    kinfeats{i} = kin(i).dat;

    % motion energy
    temp = obj(i).me.data; %cell(trials,1)
    taxis = dfparams.time;
    if dfparams.warp
        alignTimes = obj(i).bp.ev.([dfparams.alignEv '_warp']);
    else
        alignTimes = obj(i).bp.ev.(dfparams.alignEv);
    end
    tempme = zeros(numel(taxis),numel(temp));
    for j = 1:numel(temp)
        if ~dfparams.warp
            tempme(:,j) = interp1(obj(i).traj{1}(j).frameTimes-0.5-alignTimes(j),temp{j},taxis);
        else
            tempme(:,j) = interp1(obj(i).traj{1}(j).frameTimes_warped-alignTimes(j),temp{j},taxis);
        end
    end

    kinfeats{i} = cat(3, kinfeats{i}, tempme);
    kin(i).featLeg{end+1} = 'motion_energy';


    % %     STANDARDIZE FEATURES (ZERO MEAN, UNIT VARIANCE)
    kinfeats{i} = standardizeFeatures(kinfeats{i});
    kinfeats{i}(1:5,:,:) = kinfeats{i}(1:5,:,:) * 0;



end

end

