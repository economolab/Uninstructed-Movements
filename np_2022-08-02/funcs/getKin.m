function [kin,kinfeats,kinfeats_reduced] = getKin(meta,obj,dfparams,params,me)

% getKinematicsFromVideo() returns 2 variables
% - kin: feature matrix of size (time,trials,features). Features defined in
%         params.traj_features
% - featLeg: legend corresponding to features in kin (for 2nd dimension)
for i = 1:numel(meta)
    [kin(i).dat,kin(i).featLeg] = getKinematicsFromVideo(obj(i),dfparams,params(i),1:obj(i).bp.Ntrials);
    
%     % TONGUE ANGLE AND LENGTH % TODO
    [ang,len] = getLickAngleAndLength(kin(i).featLeg,kin(i).dat);
    kin(i).featLeg{end+1} = 'tongue_angle';
    kin(i).featLeg{end+1} = 'tongue_length';
    
    kinfeats{i} = kin(i).dat;
%     % create feature matrix, feats, and assign to a field in dat
    kinfeats{i} = cat(3,kin(i).dat,reshape(ang,size(ang,1),size(ang,2),1));
    kinfeats{i} = cat(3,kinfeats{i},reshape(len,size(len,1),size(len,2),1));
    
    
    % MOTION ENERGY
    if me(i).use
        kinfeats{i} = cat(3,kinfeats{i},reshape(me(i).data,size(me(i).data,1),size(me(i).data,2),1));
        kin(i).featLeg{end+1} = 'motion_energy';
    end
    % To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m
    
    
    % STANDARDIZE FEATURES (ZERO MEAN, UNIT VARIANCE)
    kinfeats{i} = standardizeFeatures(kinfeats{i});
    kinfeats{i}(1:5,:,:) = kinfeats{i}(1:5,:,:) * 0;
    
    
    % DIMENSIONALITY REDUCTION
    % many of the video features will be highly correlated, so we will perform PCA/FA
    % on the matrix of features to reduce the dimensionality to a set of factors that
    % best explain the movement captured by the video recordings
    kinfeats_reduced{i} = reduceDimensionVideoFeatures(kinfeats{i},dfparams.feat_varToExplain);
    
end



disp('DONE CREATING FEATURE MATRIX AND REDUCED DIM FEATURE MATRIX')

end

