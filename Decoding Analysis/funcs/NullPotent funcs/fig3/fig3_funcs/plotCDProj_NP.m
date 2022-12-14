function plotCDProj_NP(allrez_potent,allrez_null,rez_potent,rez_null,sav,titlestring,plotmiss)
pptx.newVersions = [1 0 0]; % 1 - creates new version of pptx, 0 - add to existing if already exists, 1 entry for each figure created here


clrs = getColors();
lw = 1;
alph = 0.12;

sample = mode(rez_potent(1).ev.sample - rez_potent(1).align);
delay = mode(rez_potent(1).ev.delay - rez_potent(1).align);


for i = 1:numel(rez_potent(1).cd_labels) % for each coding direction
    f = figure; hold on
%     f.Position = [680   748   396   230];
    ax = gca;
%         ax = nexttile; hold on;
    tempdat = squeeze(allrez_potent.cd_proj(:,:,i,:));
    tempmean = mean(tempdat,3);
    temperror = std(tempdat,[],3)./numel(rez_potent);
    shadedErrorBar(rez_potent(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(rez_potent(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    if plotmiss
        sm = 21;
        shadedErrorBar(rez_potent(1).time,mySmooth(tempmean(:,3),sm),mySmooth(temperror(:,3),sm),{'Color',clrs.rhit*0.5,'LineWidth',lw},alph, ax)
        shadedErrorBar(rez_potent(1).time,mySmooth(tempmean(:,4),sm),mySmooth(temperror(:,4),sm),{'Color',clrs.lhit*0.5,'LineWidth',lw},alph, ax)
    end
    tempdat = squeeze(allrez_null.cd_proj(:,:,i,:));
    tempmean = mean(tempdat,3);
    temperror = std(tempdat,[],3)./sqrt(numel(rez_null));
    shadedErrorBar(rez_potent(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit_aw,'LineWidth',lw,'LineStyle','-'},alph, ax)
    shadedErrorBar(rez_potent(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit_aw,'LineWidth',lw,'LineStyle','-'},alph, ax)
    if plotmiss
        sm = 21;
        shadedErrorBar(rez_null(1).time,mySmooth(tempmean(:,3),sm),mySmooth(temperror(:,3),sm),{'Color',clrs.rhit_aw*0.5,'LineWidth',lw},alph, ax)
        shadedErrorBar(rez_null(1).time,mySmooth(tempmean(:,4),sm),mySmooth(temperror(:,4),sm),{'Color',clrs.lhit_aw*0.5,'LineWidth',lw},alph, ax)
    end

%     xlim([rez(1).time(1);rez(1).time(end)])
    xlim([rez_potent(5).time(1);2])

    title([rez_null(1).cd_labels{i} ' | ' titlestring])
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax.FontSize = 12;

    xline(sample,'k:','LineWidth',2)
    xline(delay,'k:','LineWidth',2)
    xline(0,'k:','LineWidth',2)

    curmodename = rez_null(1).cd_labels{i};
    shadetimes = rez_null(1).time(rez_null(1).cd_times.(curmodename));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    %     y = [-60 -60 50 50];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';


    ylim([y(1) y(3)]);
    %     ylim([-60 50]);

    if sav
%         pptx.figPath = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\figs\st_elsayed\cd'; % path to save pptx file to
%         pptx.filename = ['st_elsayed_' spacename '_cd']; % name of ppt file, if already exists, will add fig to new slide
%         pptx.newVersion = pptx.newVersions(i);
%         pptx.slideTitle = [rez(1).cd_labels{i} ' | ' spacename];
%         pptx.fig = f;
%         myExportToPPTX(pptx)

        pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\fig1_v2\figs\cd1';
        fn = ['cd_' lower(fns{i}(3:end-7))];
        mysavefig(f,pth,fn);
%                 exportfig(f, fullfile(pth,fn),'Format','eps','Color','rgb')
    end

    hold off

end


end