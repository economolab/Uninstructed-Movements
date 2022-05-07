function plotAllCDs(rez,ev,alignEv,plt,toPlot)

if nargin > 4
    rez.psth = rez.(toPlot)
end

% get field names for each mode
[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');

psth = rez.psth;
datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
datacov(isnan(datacov)) = 0;
eigsum = sum(eig(datacov));

if strcmpi(alignEv,'gocue')
    sample = mode(ev.sample) - mode(ev.(alignEv));
    delay  = mode(ev.delay) - mode(ev.(alignEv));
else
    % if there are lots of no response trials, mode=0
    % get estimate of best alignment time with median
    sample = mode(ev.sample) - median(ev.(alignEv));
    delay  = mode(ev.delay) - median(ev.(alignEv));
end
zeroEv  = 0; % corresponds to either gocue or moveonset

% plot each mode
fig = figure;
set(fig,'Position',[-1919           1        1920        1007])
for i = 1:numel(fns)
    subplot(2,1,i);
    hold on
    for j = 1:numel(plt.conditions)
        cond = plt.conditions(j);
        latent = mySmooth(psth(:,:,cond)*rez.(fns{i}),plt.smooth);
        plot(rez.time,latent,'Color',plt.colors{j},'LineWidth',plt.lw(j));
    end
    
    varExp = var_proj(rez.(fns{i}), datacov, eigsum);

    
    xline(sample,'k--','LineWidth',0.5);
    xline(delay,'k--','LineWidth',0.5);
    xline(zeroEv,'k--','LineWidth',0.5);
%     title([fns{i}(1:end-5) '  %VE = ' num2str(round(varExp*100,2))])
    title([fns{i}(1:end-5)])
    xlim([rez.time(1)+0.2,rez.time(end)])
    ax = gca;
    ax.FontSize = 20;
%     ax.YTick = [];
%     if i < 8
%         ax.XTick = [];
%     end
%     ax.Color = [237, 237, 237]./255;
    hold off
end

leg = legend(plt.legend,'Location','best');
% leg.Position = [0.60,0.066,0.28,0.13];
sgtitle(plt.title,'FontSize',30)

h=axes(fig,'visible','off');
h.XLabel.Visible='on';
h.YLabel.Visible='on';
ylabel(h,'Activity (a.u.)', 'FontSize', 30);

xlabel(h,['Time (s) from ' alignEv], 'FontSize', 30);

% matlab has an issue when using saveas to save a figure
% the axis color gets saved as white no matter what. this is the workaround
% doesn't work if you want no background color, like for an svg
set(fig,'Color',[1 1 1]); set(fig,'InvertHardCopy','off'); 

% save
if plt.save
    savepth = fullfile(pwd,plt.title);
    saveas(fig,savepth,'png')
end

end % plotAllModes





