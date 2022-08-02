function plotKinfeats(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)



[~,mask] = patternMatchCellArray(kin(1).featLeg,feats2plot,'any');
featix = find(mask);


for i = 1:numel(obj) % for each session
    align = mode(obj(i).bp.ev.(dfparams.alignEv));
    sample = mode(obj(i).bp.ev.sample) - align;
    delay = mode(obj(i).bp.ev.delay) - align;

    for k = 1:numel(featix)
        f = figure; hold on;
        f.Position = [680    85   821   893];
        t = tiledlayout('flow');

        for j = 1:numel(cond2plot)
            ax = nexttile; hold on;
            trials = params(i).trialid{cond2plot(j)};
            temp = kinfeats{i}(:,trials,featix(k));
            imagesc(dfparams.time, 1:size(temp,2), temp');
            xline(delay,'w--');
            xline(sample,'w--');
            xline(0,'w--');
            title(dfparams.cond{cond2plot(j)});
            xlim([dfparams.times(1) dfparams.times(2)]);
            ylim([0 size(temp,2)]);
        end
        xlabel(t,['time (s) from ' dfparams.alignEv])
        sgtitle([meta(i).anm ' ' meta(i).date ' | ' feats2plot{k}], 'Interpreter','none')

        if sav
           pth = [ 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\mc_stim\figs\' meta(i).anm '\kin'];
           fn = [ meta(i).anm '_' meta(i).date '_' feats2plot{k} ];
           mysavefig(f,pth,fn)
        end

    end
    
end


end