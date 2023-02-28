function plotCDProj_Ramping(allrez,rez,spacename)

if strcmp(spacename,'Null')
    clrs.Rhit = [0.5 0.5 1]; clrs.Lhit = [1 0.5 0.5];
else
    clrs.Rhit = [0 0 1]; clrs.Lhit = [1 0 0];
end

lw = 3.5;
alph = 0.5;

sample = mode(rez(1).ev.sample - rez(1).align);
delay = mode(rez(1).ev.delay - rez(1).align);


for i = 1:numel(rez(1).cd_labels) % for each coding direction
    subplot(2,2,i)
    hold on
    ax = gca;
    tempdat = squeeze(allrez.cd_proj(:,:,i,:));
    tempmean = mean(tempdat,3);
    temperror = std(tempdat,[],3)./sqrt(numel(rez));
    shadedErrorBar(rez(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.Rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(rez(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.Lhit,'LineWidth',lw},alph, ax)

%     xlim([rez(1).time(1);rez(1).time(end)])
    xlim([rez(1).time(1);2])

    title([rez(1).cd_labels{i} ' | ' spacename])
    xlabel('Time (s) from goCue')
    ylabel('Activity (a.u.)')
    ax.FontSize = 12;

    xline(sample,'k:','LineWidth',2)
    xline(delay,'k:','LineWidth',2)
    xline(0,'k:','LineWidth',2)
    %xlim([-2.5 0])
    %ylim([-0.2 0.15])

    %ylim([y(1) y(3)]);
    hold off

end


end