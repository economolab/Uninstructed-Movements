function plotMotionEnergy(obj,me,sessionix,kinfeats,kin)



sessix = sessionix;


jawpos = kinfeats{sessix}(:,:,14);
jawpos = jawpos(:);

medata = me(sessix).data;
medata = medata(:);

newtime = (1:numel(medata))./200;

ix = me(sessix).moveIx(:);


figure; 

ax(1) = subplot(2,1,1); hold on;
patchline(newtime,medata ./ 35,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
% ix = medata > me(sessix).moveThresh;
z = medata;
z(~ix) = nan;
plot(newtime,z./35,'r','LineWidth',2)

ax(2) = subplot(2,1,2); hold on;
patchline(newtime,jawpos,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
z = jawpos;
z(~ix) = nan;
plot(newtime,z,'b','LineWidth',2)

linkaxes(ax);





end






