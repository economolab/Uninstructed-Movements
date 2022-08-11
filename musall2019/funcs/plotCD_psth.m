function plotCD_psth(rez,obj,params,times,meta)

% plot projections onto coding directions 


close all
clrs = getColors();
cols{1} = clrs.rhit;
cols{2} = clrs.lhit;
lw = 6;
alph = 0.5;

sample = mode(obj(1).bp.ev.sample - obj(1).bp.ev.(params(1).alignEvent));
delay = mode(obj(1).bp.ev.delay - obj(1).bp.ev.(params(1).alignEvent));

fns = patternMatchCellArray(fieldnames(rez(1).cd.null),{'latent'},'all');




for i = 1:numel(rez) % for each session, plot null and potent cd latents
    fnull = figure;
    fnull.Position = [354   159   570   811];
    fpotent = figure;
    fpotent.Position = [927   158   595   814];

    reztemp = rez(i);
    for j = 1:numel(fns) % plot each cd
        null = reztemp.cd.null.(fns{j});
        potent = reztemp.cd.potent.(fns{j});
        
        set(0,'CurrentFigure',fnull); 
        axnull = subplot(3,1,j); hold on;

        for k = 1:size(null,2)
            plot(obj(i).time,mySmooth(null(:,k),21),'Color',cols{k},'LineWidth',2)
        end
        

        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)

        curmodename = fns{j};
        timefieldname = [lower(curmodename(3:end-7))];
        shadetimes = obj(1).time(times.(timefieldname));
        x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
        y = [axnull.YLim(1) axnull.YLim(1) axnull.YLim(2) axnull.YLim(2)];
        %     y = [-60 -60 50 50];
        fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
        fl.FaceAlpha = 0.3;
        fl.EdgeColor = 'none';
        
        sgtitle(['Null | ' meta(i).anm ' ' meta(i).date])

        set(0,'CurrentFigure',fpotent);
        axpotent = subplot(3,1,j); hold on;

        for k = 1:size(potent,2)
            plot(obj(i).time,mySmooth(potent(:,k),21),'Color',cols{k},'LineWidth',2)
        end

        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)


        curmodename = fns{j};
        timefieldname = [lower(curmodename(3:end-7))];
        shadetimes = obj(1).time(times.(timefieldname));
        x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
        y = [axpotent.YLim(1) axpotent.YLim(1) axpotent.YLim(2) axpotent.YLim(2)];
        %     y = [-60 -60 50 50];
        fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
        fl.FaceAlpha = 0.3;
        fl.EdgeColor = 'none';

        sgtitle(['Potent | ' meta(i).anm ' ' meta(i).date])


    end
end


end