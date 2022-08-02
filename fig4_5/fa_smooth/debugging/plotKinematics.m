function plotKinematics(kin,kinfeats,obj,sessionix)

six = sessionix; % session
temp = kinfeats{six}; % time,trials,feats
leg = kin(six).featLeg;

figure(1);
for i = 1:size(temp,2)
    clf
    imagesc(obj(six).time,1:numel(leg),squeeze(temp(:,i,:))');
    title(['Trial ' num2str(i)])
    set(gca,'YTick',1:numel(leg));
    set(gca,'YTickLabel',leg);
    colorbar
    pause
end



end