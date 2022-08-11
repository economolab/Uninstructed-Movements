function ix = getMoveIndex(me,kin,kinfeats)


medata = me.data;

% get jaw pos/vel view1
[~,ix] = patternMatchCellArray(kin.featLeg,{'jaw_ydisp_view1','jaw_yvel_view1'},'any');
ix = find(ix);

jawpos = kinfeats(:,:,ix(1));

ix1 = medata > me.moveThresh;
ix2 = jawpos > 0;
ix = ix1 | ix2;

% % newtime = (1:numel(medata))./200;
% % 
% % figure;
% % 
% % ax(1) = subplot(2,1,1); hold on;
% % patchline(newtime,medata ./ 35,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
% % % ix = medata > me(sessix).moveThresh;
% % z = medata;
% % z(~ix) = nan;
% % plot(newtime,z./35,'r','LineWidth',2)
% % 
% % ax(2) = subplot(2,1,2); hold on;
% % patchline(newtime,jawpos,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
% % % ix = medata > me(sessix).moveThresh;
% % z = jawpos;
% % z(~ix) = nan;
% % plot(newtime,z,'b','LineWidth',2)
% % % z = jawpos;
% % % z(z<0) = nan;
% % % plot(newtime,z,'g','LineWidth',2)
% % 
% % linkaxes(ax)


end