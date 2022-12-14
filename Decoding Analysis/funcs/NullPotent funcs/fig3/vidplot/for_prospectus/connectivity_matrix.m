
m = 256;

cmap = redblue(m);

rhsv(1,:) = [347/360 100/100 100/100];
for i = 1:99
    rhsv(i+1,:) = [rhsv(1,1) rhsv(i,2).*1 rhsv(i,3).*0.97];
end
rhsv = flip(rhsv);

bhsv = [208/360 56/100 31/100];

f = figure;
hold on;
for i = 1:size(cmap,1)
    ff = fill([0 1 1 0] , [i i i+1 i+1],cmap(i,:));
    ff.EdgeColor = 'none';
end
ylim([0 256])
ax = gca;
ax.XTick = [];
ax.YTick= [];
box off