function plotPerformanceAllMice(meta,obj,rez,dfparams,params)

perf = rez(1).perf;
for i = 2:numel(rez)
    perf = cat(1,perf,rez(i).perf);
end
perf = perf .* 100; 


f = figure; 
f.Position = [316          73        1232         905];
ax = axes(f);
violincols = reshape(cell2mat(dfparams.plt.color),3,numel(dfparams.cond))';
vs = violinplot(perf,dfparams.cond,...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols);
ylabel('Performance (%)')
ylim([0,100])
title('Performance on all sessions, all mice')
ax.FontSize = 20;


end