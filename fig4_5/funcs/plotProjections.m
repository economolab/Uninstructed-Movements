function plotProjections(params,obj,rez)

close all
clear cols clrs

sav = 0;

cols = getColors();
clrs{1} = cols.rhit;
clrs{2} = cols.lhit;
lw = 3;
alph = 0.5;
for i = 1:numel(rez)
    %     optimization_plots(rez,obj,dat,params); % old
    
    temp = rez(i).N_potent;
    
    f = figure;
    f.Position = [-1143         -27         357         848];
    for dimix = 1:size(temp,3)
        ax = subplot(size(temp,3),1,dimix); hold on
        for j = 1:2
            tempdat = squeeze(temp(:,params(i).trialid{j+1},dimix));
            means = mean(tempdat,2);
            stderr = std(tempdat,[],2) ./ numel(params(i).trialid{j+1});
            shadedErrorBar(obj(1).time,means,stderr,{'Color',clrs{j},'LineWidth',lw},alph, ax)
%             plot(obj(1).time,means,'Color',clrs{j},'LineWidth',1)
        end
        title(['Potent ' num2str(dimix)])
        xlim([obj(1).time(15),obj(1).time(end)])
        
        align = mode(obj(i).bp.ev.(params(i).alignEvent));
        sample = mode(obj(i).bp.ev.sample) - align;
        delay = mode(obj(i).bp.ev.delay) - align;
        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)
        
        ax.FontSize = 20;
        hold off;
    end
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/potent';
        fn = [meta(i).anm '_' meta(i).date];
        mysavefig(f,pth,fn);
        pause(3)
    end
    
    temp = rez(i).N_null;
    
    f = figure;
    f.Position = [-1501         187         357         631];
    for dimix = 1:size(temp,3)
        ax = subplot(size(temp,3),1,dimix); hold on
        for j = 1:2
            tempdat = squeeze(temp(:,params(i).trialid{j+1},dimix));
            means = mean(tempdat,2);
            stderr = std(tempdat,[],2) ./ numel(params(i).trialid{j+1});
            shadedErrorBar(obj(1).time,means,stderr,{'Color',clrs{j},'LineWidth',lw},alph, ax)
%             plot(obj(1).time,means,'Color',clrs{j},'LineWidth',1)
        end
        title(['Null ' num2str(dimix)])
        xlim([obj(1).time(15),obj(1).time(end)])
        
        align = mode(obj(i).bp.ev.(params(i).alignEvent));
        sample = mode(obj(i).bp.ev.sample) - align;
        delay = mode(obj(i).bp.ev.delay) - align;
        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)
        
        ax.FontSize = 20;
        hold off;
    end
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace/null';
        fn = [meta(i).anm '_' meta(i).date];
        mysavefig(f,pth,fn);
        pause(3)
    end
    
    
end