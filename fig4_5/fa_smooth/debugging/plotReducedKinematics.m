function plotReducedKinematics(kinfeats_reduced,obj,sessionix)

six = sessionix; % session
temp = kinfeats_reduced{six}; % time,trials,feats

figure(1);
for i = 1:size(temp,2)
    clf
    imagesc(obj(six).time,1:size(temp,3),squeeze(temp(:,i,:))');
    title(['Trial ' num2str(i)])
    colorbar
    pause
end



end