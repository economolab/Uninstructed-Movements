function rhos = calcCorrCoef(V,N,W)

% return correlation between each feature in V and each estimated feature
% in V_hat
V_hat = N*W;
rhos = corr(V(:),V_hat(:));

end