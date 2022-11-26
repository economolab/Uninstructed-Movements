function cd = calcCD(rez,times)
tempdat = rez.psth(:,:,[1,2]);
mu = squeeze(mean(tempdat(times,:,:),1));
sd = squeeze(std(tempdat(times,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
end