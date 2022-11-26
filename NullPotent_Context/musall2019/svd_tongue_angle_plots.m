

% svd heatmap
figure;
trix = [trials.right;trials.left];
imagesc(tm,1:nrDims,squeeze(absVidR_Cam1(:,trix,4))');
xlabel('Time (s) from go cue')
ylabel('Trial number')
clim([-4 4])
c = colorbar;
c.Label.String = "Bottom Cam SV 4";
ax = gca;
ax.FontSize = 15;

% trial-averaged svd

temp1 = mean(squeeze(vidR_Cam1(:,trials.right,4)),2);
temp2 = mean(squeeze(vidR_Cam1(:,trials.left,4)),2);

figure; hold on
plot(tm,temp1,'b','LineWidth',3);
plot(tm,temp2,'r','LineWidth',3);
xlabel('Time (s) from go cue')
ylabel('Bottom Cam SV4, Trial-Averaged')
ax = gca;
ax.FontSize = 15;

% tongue angle
angix = find(ismember(kin.featLeg,'tongue_angle'));
temp = kinfeats{1}(:,params.trials2use,angix);

figure;
trix = [trials.right;trials.left];
imagesc(tm,1:nrDims,temp(:,trix)');
xlabel('Time (s) from go cue')
ylabel('Trial number')
% clim([-4 4])
% c = colorbar;
% c.Label.String = "Tongue Angle";
ax = gca;
ax.FontSize = 15;


