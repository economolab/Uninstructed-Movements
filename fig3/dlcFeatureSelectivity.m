% first ran kinematicNullSpace and used dat.feats_reduced here for a single
% session

% this is not a good idea. We can average behavior across trials and i've
% saved one figure of that. Need to chagne code to do that again though.

temp = dat.feats_reduced; % low-d movement features on single trials

clear condition
condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
trialid = findTrials(obj,condition);

mask = ismember(dat.trials,trialid{1});
trialid{1} = find(mask);
mask = ismember(dat.trials,trialid{2});
trialid{2} = find(mask);

dataR = temp(:,trialid{1},:);
dataL = temp(:,trialid{2},:);

rtemp = reshape(dataR,size(dataR,1)*size(dataR,2),size(dataR,3));
[~,mu,sigma] = zscore(rtemp);

a = mu + 0.5*sigma;
b = mu - 0.5*sigma;
mask = (rtemp>a);
rrtemp = nan(size(rtemp));
rrtemp(mask) = 1;

t = reshape(rrtemp,size(dataR,1),size(dataR,2),size(dataR,3));
tmeantrials = squeeze(nanmean(t,2));

figure; plot(tmeantrials)


%%
k = 10;
nTrials = round(min(numel(trialid{1}),numel(trialid{2})) / 2);
corr_matrix_selectivity = zeros(size(temp,1),size(temp,1),k);
for j = 1:k
    
    dataR = temp(:,randsample(trialid{1},nTrials),:);
    dataL = temp(:,randsample(trialid{2},nTrials),:);
    dataR = reshape(dataR,size(dataR,1)*size(dataR,2),size(dataR,3));
    dataL = reshape(dataL,size(dataL,1)*size(dataL,2),size(dataL,3));
    
    selectivity = dataR - dataL;
        
    corr_matrix_selectivity = zeros(size(selectivity,1),size(selectivity,1));
    
    for i = 1:size(corr_matrix_selectivity,1)
        for j = 1:size(corr_matrix_selectivity,1)
            cc = corrcoef(selectivity(i,:),selectivity(j,:));
            corr_matrix_selectivity(i,j,k) = cc(1,2);
        end
    end
    
    disp(i)
    
end

corr_matrix_selectivity = mean(corr_matrix_selectivity,3);

%%


sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay  = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

figure; hold on;
imagesc(obj.time,obj.time,corr_matrix_selectivity);
colorbar; caxis([0 max(max(corr_matrix_selectivity))]);

lw = 4;
xline(sample,'w--','LineWidth',lw); yline(sample,'w--','LineWidth',lw)
xline(delay,'w--','LineWidth',lw); yline(delay,'w--','LineWidth',lw)
xline(0,'w--','LineWidth',lw); yline(0,'w--','LineWidth',lw)

xlim([obj.time(1)+0.2,obj.time(end)]);
ylim([obj.time(1)+0.2,obj.time(end)])

ax = gca;
ax.FontSize = 20;
hold off
colormap(hot)


