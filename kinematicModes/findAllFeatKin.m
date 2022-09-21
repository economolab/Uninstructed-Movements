function kin = findAllFeatKin(params, taxis, obj,conditions,met)
if strcmp(params.kinfind,'pos')
    view = 2; % side
    feat = 8; % jaw
    [kin.jawPos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
    view = 1; % side
    feat = 7; % nose
    [kin.nosePos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
    view = 2; % side
    feat = 1; % tongue
    [kin.tonguePos,~] = getTongueKinematics(taxis,obj,conditions,met,view,feat,params);
%     view = 2; % bottom
%     feat = 5; % top_paw
%     [kin.topPawPos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
%     view = 2; % bottom
%     feat = 6; % bottom_paw
%     [kin.bottomPawPos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
elseif strcmp(params.kinfind,'vel')
    view = 1; % side
    feat = 4; % jaw
    [~,kin.jawVel] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
    view = 1; % side
    feat = 7; % nose
    [~,kin.noseVel] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
    view = 1; % bottom
    feat = 1; % top_tongue
    [~,kin.tongueVel] = getTongueKinematics(taxis,obj,conditions,met,view,feat,params);
    %         view = 2; % bottom
    %         feat = 5; % top_paw
    %         [~,kin.topPawVel] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
    %         view = 2; % bottom
    %         feat = 6; % bottom_paw
    %         [~,kin.bottomPawVel] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params);
end
end