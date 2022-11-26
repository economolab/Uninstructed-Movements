function plotCorrelation_Cddelay_ME_v2(obj,params,rez,me,meta,cond2use)

clrs = getColors();
cols{1} = clrs.rhit;
cols{2} = clrs.lhit;

r2 = [];
% avgcd = [];
% avgme = [];


for sessix = 1:numel(obj)
%     if sessix ~=10
%         continue
%     end
    f = figure; hold on


    for cix = 1:numel(cond2use)
        % only use 400 ms up to go cue for correlations
        [~,ix1] = min(abs(obj(sessix).time - -0.4));
        [~,ix2] = min(abs(obj(sessix).time - 0));
        trix = params(sessix).trialid{cond2use(cix)};
        cddelayix = find(ismember( rez(sessix).cd_labels,'late'));
        plotcd = rez(sessix).cd_proj_trialdat(ix1:ix2,trix,cddelayix);
        avgcd =  mean(plotcd,1)';

        plotme = me(sessix).data(ix1:ix2,trix);
        avgme = mean(plotme,1)';

        mdl{cix} = fitlm(avgme,avgcd);
        % mdl.Coefficients.Estimate
        %  y ~ 1 + x1
        ax = plot(mdl{cix}); hold on;
        scatter(avgme,avgcd,70,'MarkerEdgeColor','w','MarkerFaceColor',cols{cix},'MarkerFaceAlpha',1); hold on;
        ax(1).XData = []; % remove x's plotted by ax=plot(mdl)
        ax(1).YData = [];
        ax(1) = [];
        ax(1).LineWidth = 0.5;
        ax(1).Color = cols{cix};
        ax(2).Color = cols{cix};
        ax(3).Color = cols{cix};
        leg = legend();
        set(leg,'Visible','off')

    end
    title([meta(sessix).anm '_' meta(sessix).date '  |  $R1^2$ = ' num2str(mdl{1}.Rsquared.Ordinary) '  |   $R2^2$ = ' num2str(mdl{2}.Rsquared.Ordinary)],'Interpreter','latex')
    xlabel('Late Delay Motion Energy')
    ylabel('Late Delay Coding Direction')
end



% for sessix = 1:numel(rez)
%
%     % only use 400 ms up to go cue for correlations
%     [~,ix1] = min(abs(obj(sessix).time - -0.4));
%     [~,ix2] = min(abs(obj(sessix).time - 0));
%
%     for cix = 1:numel(cond2use)
%         trix = params(sessix).trialid{cond2use(cix)};
%         cddelayix = find(ismember( rez(sessix).cd_labels,'late'));
%         plotcd = rez(sessix).cd_proj_trialdat(ix1:ix2,trix,cddelayix);
%         avgcd = mean(plotcd,1)';
%
%         plotme = me(sessix).data(ix1:ix2,trix);
%         avgme = mean(plotme,1)';
%                 scatter(avgme,avgcd,20,'MarkerEdgeColor','w','MarkerFaceColor',cols{cix},'MarkerFaceAlpha',1)
% %         scatter(avgme,avgcd,20,'MarkerEdgeColor','w','MarkerFaceColor','k','MarkerFaceAlpha',1)
%         hold on
%     end
%
%
% end
%
% f = figure; hold on
% scatter(avgme,avgcd,30,'k')
% histogram(r2,10,'EdgeColor',[29, 50, 84]./255,'FaceColor',[62, 109, 184]./255, 'FaceAlpha',0.6)



end