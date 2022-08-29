clear,clc,close all

%%
pth = 'C:\Users\munib\Documents\Economo-Lab\data\';

% load neural data and params
fn = 'JEB7_2021-04-29_run4';
% fn = 'JGR2_2021-11-17_run4';
load(fullfile(pth,'fa',fn))

% load kinematics
fn = 'kin_JEB7_2021-04-29_run4';
% fn = 'kin_JGR2_2021-11-17_run4';
load(fullfile(pth,'kinematics',fn))


%%  motion energy

trix = [params.trialid{2} ; params.trialid{3}]; % right;left trials

ix = find(ismember(kin.featLeg,{'motion_energy'}));

rez.me = kinfeats(:,trix,ix);
rez.me_tavg(:,1) = mean(kinfeats(:,params.trialid{2},ix),2);
rez.me_tavg(:,2) = mean(kinfeats(:,params.trialid{3},ix),2);

%%
figure;
subplot(1,2,1)
imagesc(obj.time,1:numel(trix),rez.me')
subplot(1,2,2); hold on
plot(obj.time,rez.me_tavg(:,1),'b','LineWidth',2)
plot(obj.time,rez.me_tavg(:,2),'r','LineWidth',2)
xlabel('Time (s) from go cue')
ylabel('Motion Energy')
xline(0,'k--','LineWidth',2)
ax = gca;
ax.FontSize = 14;

%% find a motor planning CD


e1 = mode(obj.bp.ev.goCue) - 0.5 - mode(obj.bp.ev.(params.alignEvent));
e2 = mode(obj.bp.ev.goCue) - 0.1 - mode(obj.bp.ev.(params.alignEvent));

times.late = (obj.time>e1 & obj.time<e2)';
rez.cdlate = calcCD(obj.psth(:,:,[2 3]),times.late);

% trial-averaged projections
rez.cdlate_psth(:,1) = obj.psth(:,:,1) * rez.cdlate; 
rez.cdlate_psth(:,2) = obj.psth(:,:,2) * rez.cdlate; 

% single trial projections
for i = 1:size(obj.trialdat,3)
    rez.cdlate_trialdat(:,i) = obj.trialdat(:,:,i) * rez.cdlate;
end
rez.cdlate_trialdat = rez.cdlate_trialdat(:,trix);


figure;
subplot(1,2,1)
imagesc(rez.cdlate_trialdat')
subplot(1,2,2)
plot(rez.cdlate_psth)

%% plot avg motion energy in late delay against avg cdlate projection in late delay

me_latedelay = mean(rez.me(times.late,:),1);

cd_latedelay = mean(rez.cdlate_trialdat(times.late,:),1);

mdl = fitlm(me_latedelay,cd_latedelay);
% mdl.Coefficients.Estimate
%  y ~ 1 + x1

f = figure; hold on;
ax = plot(mdl);
scatter(me_latedelay(1:numel(params.trialid{2})) , cd_latedelay(1:numel(params.trialid{2})) , 'b','filled')
scatter(me_latedelay((numel(params.trialid{2})+1:end)) , cd_latedelay((numel(params.trialid{2})+1:end)) , 'r','filled')
ax(2).Color = 'k';
ax(3).Color = 'k';
ax(4).Color = 'k';
leg = legend();
set(leg,'Visible','off')
title(['$R^2$ = ' num2str(mdl.Rsquared.Ordinary) ' $\|$ p = ' num2str(mdl.Coefficients.pValue(2))],'Interpreter','latex')
xlabel('Late Delay Motion Energy')
ylabel('Late Delay Coding Direction')
ax = gca;
ax.FontSize = 14;



%% cd at every time point

% cd = (t,n)
% trialdat = (t,n,trials)


times.all = logical(eye(numel(obj.time)));
for i = 1:numel(obj.time)
    rez.cdall(i,:) = calcCD_everyTimePoint(obj.psth(:,:,[2 3]),times.all(i,:));

    rez.cdall_trialdat(i,:) =  rez.cdall(i,:) * squeeze(obj.trialdat(i,:,:));
end
rez.cdall_psth(:,1) = mean(rez.cdall_trialdat(:,params.trialid{2}),2);
rez.cdall_psth(:,2) = mean(rez.cdall_trialdat(:,params.trialid{3}),2);
rez.cdall_trialdat = rez.cdall_trialdat(:,trix);


figure;
subplot(1,2,1)
imagesc(rez.cdall_trialdat')
subplot(1,2,2)
plot(rez.cdall_psth)

%% plot avg motion energy in late delay against avg cd projection in late delay

me_latedelay = mean(rez.me(times.late,:),1);

cd_latedelay = mean(rez.cdall_trialdat(times.late,:),1);

mdl = fitlm(me_latedelay,cd_latedelay);
% mdl.Coefficients.Estimate
%  y ~ 1 + x1

f = figure; hold on;
ax = plot(mdl);
scatter(me_latedelay(1:numel(params.trialid{2})) , cd_latedelay(1:numel(params.trialid{2})) , 'b','filled')
scatter(me_latedelay((numel(params.trialid{2})+1:end)) , cd_latedelay((numel(params.trialid{2})+1:end)) , 'r','filled')
ax(2).Color = 'k';
ax(3).Color = 'k';
ax(4).Color = 'k';
leg = legend();
set(leg,'Visible','off')
title(['$R^2$ = ' num2str(mdl.Rsquared.Ordinary) ' $\|$ p = ' num2str(mdl.Coefficients.pValue(2))],'Interpreter','latex')
xlabel('Late Delay Motion Energy')
ylabel('Late Delay Coding Direction')
ax = gca;
ax.FontSize = 14;





%% Helper functions


function cd = calcCD(dat,times)

mu = squeeze(mean(dat(times,:,:),1));

sd = squeeze(std(dat(times,:,:),[],1));

cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)

end

function cd = calcCD_everyTimePoint(dat,times)

mu = squeeze(mean(dat(times,:,:),1));


cd = ((mu(:,1)-mu(:,2)));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)

end















