function plotChoiceMode_HazardedDel(params,latentChoice,colors,taxis)
for g=1:length(params.delay)
    subplot(2,3,g)
    plot(taxis,latentChoice.left{g},'Color',colors{2},'LineWidth',2)
    hold on;
    plot(taxis,latentChoice.right{g},'Color',colors{1},'LineWidth',2)
    xlabel('Time since delay onset')
    ylabel('a.u.')
    xlim([-2 2])
    xline(0,'LineStyle','--','LineWidth',1.3)
    xline(-1.3,'LineStyle',':','LineWidth',1.3)
    xline(params.delay(g),'LineStyle','-.','LineWidth',1.3)
    legend('Left','Right','Delay onset','Sample onset','GoCue','Location','best')
    len = num2str(params.delay(g));
    subtitle=strcat('Delay length =',{' '},len);
    title(subtitle)
end
end