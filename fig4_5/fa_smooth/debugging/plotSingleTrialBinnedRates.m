function plotSingleTrialBinnedRates(obj,sessionix)

six = sessionix; % session
temp = obj(six).trialdat; % time,feats,trials

figure(1);
for i = 1:size(temp,3)
    clf
    imagesc(obj(six).time,1:size(temp,2),squeeze(temp(:,:,i))');
%     set(gca,'YTick',1:size(temp,2));
    xlabel('time (s)')
    ylabel('neurons')
    title(['Trial ' num2str(i)]);
    xline(0,'w')
    colorbar;
    pause
end

end