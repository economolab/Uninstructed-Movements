function plot_nonGC_moveTransitions_singleTrials_v5(dat,obj,me,rez,params,ndims,meta,sav)

rng(pi)


% q2m
for sessix = 1:numel(me)

    dt = params(sessix).dt;

    % transition data
    qdat = dat.qdat{sessix};
    nTrials = numel(qdat);


    % sort null/potent spaces by ve
    [~,nullix] = sort(rez(sessix).ve.null,'descend');
    nullix = nullix(1:ndims);

    [~,potentix] = sort(rez(sessix).ve.potent,'descend');
    potentix = potentix(1:ndims);

    nbouts = 5; % approximately max number of bouts in a trial
    allnull = nan(500,nTrials*nbouts,ndims);
    allpotent = nan(500,nTrials*nbouts,ndims);
    allme = nan(500,nTrials*nbouts,ndims);

    for dimix = 1:ndims

        for trix = 1:nTrials
            qdat_trix = qdat{trix};
            if isempty(qdat_trix)
                continue
            end

            sm = 41;
            null = mySmooth(rez(sessix).N_null(:,trix,nullix(dimix)),sm);
            potent = mySmooth(rez(sessix).N_potent(:,trix,potentix(dimix)),sm);
            move = mySmooth(me(sessix).move(:,trix),sm);
            medat = me(sessix).data(:,trix);

            % for each transition
            for i = 1:size(qdat_trix,2)
                tix = qdat_trix(1,i):qdat_trix(end,i);
                ix = find(tix == 0);
                if ~isempty(ix)
                    tix(ix) = [];
                end
                ts = tix .* dt;
                ts = ts - mean(ts);

                nullplot = null(tix);
                potentplot = potent(tix);
                meplot = medat(tix);

                % get center of tix
                center = ceil(numel(tix)/2);
                try
                    allnull(250-center:250+center-1,(trix+i-1),dimix) = nullplot;
                    allpotent(250-center:250+center-1,(trix+i)-1,dimix) = potentplot;
                    allme(250-center:250+center-1,(trix+i)-1,dimix) = meplot;
                catch
                    allnull(250-center:250+center-2,(trix+i-1),dimix) = nullplot;
                    allpotent(250-center:250+center-2,(trix+i)-1,dimix) = potentplot;
                    allme(250-center:250+center-2,(trix+i)-1,dimix) = meplot;
                end


            end

        end


    end



    %% sumsqmag
    meannull = sum(allnull.^2,3,'omitnan');
    meanpotent = sum(allpotent.^2,3,'omitnan');

    meannull = normalize(meannull,'range',[0 1]);
    meanpotent = normalize(meanpotent,'range',[0 1]);

    % var
    %     meannull = squeeze(nanmean(allnull.^2,3));
    %     meanpotent = squeeze(nanmean(allpotent.^2,3));

    mask = meannull == 0;
    %     mask = isnan(meannull);
    % find cols of mask with all 0s
    temp = sum(~mask,1);
    meannull(:,temp==0) = [];

    mask = meanpotent == 0;
    %     mask = isnan(meanpotent);
    % find cols of mask with all 0s
    temp = sum(~mask,1);
    meanpotent(:,temp==0) = [];


    cm = parula;
    f = figure;
    f.Position = [668    42   293   954];

    xlims = [-0.4,0.4];

    %% me
    ax(1) = nexttile; hold on;
    toplot = allme(:,:,1);
    mask = isnan(toplot);
    % find cols of mask with all 0s
    temp = sum(~mask,1);
    toplot(:,temp==0) = [];
    toplot = normalize(toplot,'range',[0 1]);
    %
    %     ttt = toplot(250-31:250+31-2,:);
    %     figure; imagesc(ttt')

    %     mask = toplot > me(sessix).moveThresh;
    %     toplot(~mask) = -max(max(toplot));

    [~,ix] = sort(nanmean(toplot,1),'descend');

    imagesc(obj(sessix).time,1:size(toplot,2),toplot(:,ix)')
    xlim(xlims)
    ylim([0 size(toplot,2)])
    xlabel('Time to move (s)')
    ylabel('Motion Energy')
    ax(1).FontSize = 12;
    colorbar(ax(1)); %caxis([0 30])
    xline(0,'k:','LineWidth',2)
    colormap(cm)
    %% null and potent

    ax(end+1) = nexttile; hold on;
    toplot = meanpotent;
    imagesc(obj(sessix).time,1:size(toplot,2),toplot(:,ix)')
    xlim(xlims)
    ylim([0 size(toplot,2)])
    xlabel('Time to move (s)')
    ylabel('Potent')
    ax(2).FontSize = 12;
    colorbar(ax(2));% caxis([0 2.5])
    xline(0,'k:','LineWidth',2)
    colormap(cm)

    ax(end+1) = nexttile; hold on;
    toplot = meannull;
    imagesc(obj(sessix).time,1:size(toplot,2),toplot(:,ix)')
    xlim(xlims)
    ylim([0 size(toplot,2)])
    xlabel('Time to move (s)')
    ylabel('Null')
    ax(3).FontSize = 12;
    colorbar(ax(3)); %caxis([0 2.5])
    xline(0,'k:','LineWidth',2)
    colormap(cm)



    %%

    sgtitle([meta(sessix).anm ' - ' meta(sessix).date])

    if sav
        ppt.fig = f;
        ppt.figPath = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\figs\move-transitions';
        ppt.filename = 'moveToQuietSingleTrials';
        ppt.newVersion = 0;
        ppt.slideTitle = [meta(sessix).anm '  ' meta(sessix).date];
        myExportToPPTX(ppt)
    end

    clear ax




end



end