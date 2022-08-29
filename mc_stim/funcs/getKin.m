function [kin,kinfeats] = getKin(meta,obj,dfparams,params)

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
    
    
    
    % STANDARDIZE FEATURES (ZERO MEAN, UNIT VARIANCE)
    kinfeats{i} = standardizeFeatures(kinfeats{i});
    kinfeats{i}(1:5,:,:) = kinfeats{i}(1:5,:,:) * 0;
    

    
end



disp('DONE CREATING FEATURE MATRIX AND REDUCED DIM FEATURE MATRIX')

end

