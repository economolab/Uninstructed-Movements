
%% overview of data

% side cam: view = 1
% relevant coordinates: (x,z)
% % view 1 bodyparts:
%     1{'tongue'  } 
%     2{'jaw'     } 
%     3{'nose'    } 
%     4{'lickport'} 

% bottom cam: view = 2
% relevant coordinates: (x,y)
% % view 2 bodyparts:
%     1{'top_tongue'       } 
%     2{'topleft_tongue'   } 
%     3{'bottom_tongue'    } 
%     4{'leftbottom_tongue'} 
%     5{'top_paw'          } 
%     6{'bottom_paw'       } 
%     7{'lickport'         } 
%     8{'jaw'              } 

% obj.traj = 1x2 cell array
% obj.traj{1 (sidecam)} = {1xnumTrials} struct
% obj.traj{2 (bottomcam)} = {1xnumTrials} struct
%     struct fields: fn, featNames, ts
%     ts is trajectories
% obj.traj{1}(1) is struct for trial 1
% obj.traj{1}(1).ts is (time,coords,bodypart)
% sidecam has x,z coords
% bottomcam has x,y coords

