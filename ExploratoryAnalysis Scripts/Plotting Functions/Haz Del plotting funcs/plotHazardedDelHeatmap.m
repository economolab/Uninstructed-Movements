function plotHazardedDelHeatmap (met,vel_by_cond, delaylen,taxis,ix1,ix2)

% Make a heatmap of single-trial jaw velocities, sorted by delay length
figure();
subplot(1,2,2)
imagesc(taxis,1:length(met.trialid{1}),vel_by_cond{1}(:,:)')
title('Right trials')
xlabel('Time since go-cue (s)')
ylabel('Trials')
colorbar
caxis([0 10])
for ii = 1:length(met.del_trialid{1})
    tri = met.trialid{1}(ix1);
    delstart = 0-delaylen(tri(ii));
    sampstart = delstart-1.3;
    line([delstart,delstart],[ii-0.5,ii+0.5],'Color','Green')
    line([sampstart,sampstart],[ii-0.5,ii+0.5],'Color','Magenta')
end
line([0 0],[1 length(met.trialid{1})],'Color','White','LineStyle','--')

subplot(1,2,1)
imagesc(taxis,1:length(met.del_trialid{2}),vel_by_cond{2}(:,:)')
title('Left trials')
colorbar
caxis([0 10])
for ii = 1:length(met.del_trialid{2})
    tri = met.trialid{2}(ix2);
    delstart = 0-delaylen(tri(ii));
    sampstart = delstart-1.3;
    line([delstart,delstart],[ii-0.5,ii+0.5],'Color','Green')
    line([sampstart,sampstart],[ii-0.5,ii+0.5],'Color','Magenta')
end
line([0 0],[1 length(met.trialid{2})],'Color','White','LineStyle','--')

xlabel('Time since go-cue (s)')
ylabel('Trials')
end %plotHazardedDelHeatmap