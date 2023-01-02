function plotNP_CD_Ramp(allrez_null,rez_null,allrez_potent,rez_potent)

clrs.Rpotent = [0 0 1]; clrs.Lpotent = [1 0 0]; clrs.Rnull = [0.5 0.5 1]; clrs.Lnull = [1 0.5 0.5];
lw = 3.5;
alph = 0.5;

sample = mode(rez_null(1).ev.sample - rez_null(1).align);
delay = mode(rez_null(1).ev.delay - rez_null(1).align);

modeix = find(strcmp(rez_null(1).cd_labels,'ramping'));
figure();
hold on
ax = gca;
tempdat = squeeze(allrez_null.cd_proj(:,:,modeix,:));

start = find(rez_null(1).time>-2.5,1,'first');
stop = find(rez_null(1).time<sample,1,'last');
tempmean = mean(tempdat,3);
presampMu = mean(tempmean(start:stop,:),1);
% Presample mean-center
tempmean = tempmean-presampMu;
temperror = std(tempdat,[],3)./sqrt(numel(rez_null));
shadedErrorBar(rez_null(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.Rnull,'LineWidth',lw},alph, ax)
shadedErrorBar(rez_null(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.Lnull,'LineWidth',lw},alph, ax)

tempdat = squeeze(allrez_potent.cd_proj(:,:,modeix,:));
tempmean = mean(tempdat,3);
% Presample mean-center
presampMu = mean(tempmean(start:stop,:),1);
tempmean = tempmean-presampMu;
temperror = std(tempdat,[],3)./sqrt(numel(rez_potent));
shadedErrorBar(rez_potent(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.Rpotent,'LineWidth',lw},alph, ax)
shadedErrorBar(rez_potent(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.Lpotent,'LineWidth',lw},alph, ax)

%     xlim([rez(1).time(1);rez(1).time(end)])
xlim([rez_null(1).time(1);2])

xlabel('Time (s) from goCue')
ylabel('Activity (a.u.)')
ax.FontSize = 12;
title('Ramping mode | Null and Potent')

xline(sample,'k:','LineWidth',2)
xline(delay,'k:','LineWidth',2)
xline(0,'k:','LineWidth',2)


hold off

end
