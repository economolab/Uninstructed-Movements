function [ang,len] = getLickAngleAndLength(dat,kin)

ix1 = find(~cellfun(@isempty, strfind(dat.featLeg, 'top_tongue_xdisp')));
ix2 = find(~cellfun(@isempty, strfind(dat.featLeg, 'topleft_tongue_xdisp_view2')));
t1x = kin(:, :, ix1);
t2x = kin(:, :, ix2);

ix1 = find(~cellfun(@isempty, strfind(dat.featLeg, 'top_tongue_ydisp')));
ix2 = find(~cellfun(@isempty, strfind(dat.featLeg, 'topleft_tongue_ydisp_view2')));
t1y = kin(:, :, ix1);
t2y = kin(:, :, ix2);

ix1 = find(~cellfun(@isempty, strfind(dat.featLeg, 'bottom_tongue_xdisp')));
ix2 = find(~cellfun(@isempty, strfind(dat.featLeg, 'bottomleft_tongue_xdisp')));
t3x = kin(:, :, ix1);
t4x = kin(:, :, ix2);

ix1 = find(~cellfun(@isempty, strfind(dat.featLeg, 'bottom_tongue_ydisp')));
ix2 = find(~cellfun(@isempty, strfind(dat.featLeg, 'bottomleft_tongue_ydisp')));
t3y = kin(:, :, ix1);
t4y = kin(:, :, ix2);

ix1 = (~cellfun(@isempty, strfind(dat.featLeg, 'jaw_xdisp_view2')));
jx = kin(:, :, ix1);
ix1 = (~cellfun(@isempty, strfind(dat.featLeg, 'jaw_ydisp_view2')));
jy = kin(:, :, ix1);



tongue_tip_x = 0.5.*(t1x+t3x);
tongue_tip_y = 0.5.*(t1y+t3y);

tongue_base_x = 0.5.*(t2x+t4x);
tongue_base_y = 0.5.*(t2y+t4y);

dx = (tongue_tip_x-jx);
dy = (tongue_tip_y-jy);
ang = atan(dy./dx);
ang(dx<0 & dy>0) = ang(dx<0 & dy>0) + pi;     %Correction for the quadrant that the angle lies in
ang(dx<0 & dy<0) = ang(dx<0 & dy<0) - pi;
len = sqrt(dx.^2 + dy.^2);


end