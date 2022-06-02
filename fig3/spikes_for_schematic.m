f = figure(1); 
f.Position = [-1411         175         207         461];
for i = 56
    clf
    plot(obj.clu{1}(i).trialtm(1:20:end),obj.clu{1}(i).trial(1:20:end),'k.','MarkerSize',1)
    xlim([0 5])
    ylim([0 350])
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
end

pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/schematic';
fn = 'JEB7_clunum56';
mysavefig(f,pth,fn)