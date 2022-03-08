function calcEigenpose(obj)

feats = [2 3];

trix = 1;
ts = obj.traj{1}(trix).ts(:,[1 2],feats);
ts = reshape(ts,size(ts,1),4);

for i = 1:size(ts,2)
    ts(:,i) = fillmissing(ts(:,i),'nearest');
    ts(:,i) = medfilt1(ts(:,i),3);
end

[pcs, latents, eigval, ~, varexp] = pca(ts-nanmean(ts));

tm = (1:size(ts,1))./400;

figure; 
subplot(121)
plot(latents)
subplot(122)
scatter3(latents(:,1),latents(:,2),latents(:,3),20,tm); colorbar;

end