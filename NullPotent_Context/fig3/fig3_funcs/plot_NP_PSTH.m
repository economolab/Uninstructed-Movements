function plot_NP_PSTH(rez,obj,params,ndims,cond2plot,meta)

lw = 2;
clrs = getColors_Context;
cols{1} = clrs.AFChit;
cols{2} = clrs.AWhit;

align = mode(obj(1).bp.ev.(params(1).alignEvent));

sample = mode(obj(1).bp.ev.sample - align);
delay = mode(obj(1).bp.ev.delay - align);


% % - null
sm = 31;
for sessix = 1:numel(rez)
    f = figure;
    %f.Position = [92   144   610   799];

    % sort dimensions by variance explined (most to least)
    [~,ix] = sort(rez(sessix).ve.null,'descend');

    for dimix = 1:ndims
        ax = nexttile;
        hold on;
        for condix = 1:numel(cond2plot)
            plot(obj(1).time,mySmooth(rez(sessix).N_null_psth(:,ix(dimix),cond2plot(condix)),sm),"Color",cols{condix},'LineWidth',lw)
        end
        title([num2str(round(rez(sessix).ve.null(ix(dimix))*100),3) ' %VE'])
        xline(sample,'k:')
        xline(delay,'k:')
        xline(0,'k:')
        xlim([-2.6 2.5])
    end
    sgtitle(['Null | ' meta(sessix).anm ' - ' meta(sessix).date])
end

% % - potent
sm = 21;
for sessix = 1:numel(rez)
    f = figure;
    %f.Position = [917   141   610   799];

    % sort dimensions by variance explined (most to least)
    [~,ix] = sort(rez(sessix).ve.potent,'descend');

    for dimix = 1:ndims
        ax = nexttile;
        hold on;
        for condix = 1:numel(cond2plot)
            plot(obj(1).time,mySmooth(rez(sessix).N_potent_psth(:,ix(dimix),cond2plot(condix)),sm),"Color",cols{condix},'LineWidth',lw)
        end
        title([num2str(round(rez(sessix).ve.potent(ix(dimix))*100),3) ' %VE'])
        xline(sample,'k:')
        xline(delay,'k:')
        xline(0,'k:')
        xlim([-2.6 2.5])
    end
    sgtitle(['Potent | ' meta(sessix).anm ' - ' meta(sessix).date])

end




end