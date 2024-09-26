function catpsth = concatPSTH_byCond(psth)
% inputs
%  - psth: should be psth matrix of size (time,clusters,conditions)

catpsth = psth(:,:,1);
for i = 2:size(psth,3)
     catpsth = [catpsth ; psth(:,:,i)];
end