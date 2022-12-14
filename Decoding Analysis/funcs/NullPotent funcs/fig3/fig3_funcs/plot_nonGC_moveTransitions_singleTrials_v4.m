function plot_nonGC_moveTransitions_singleTrials_v4(dat,obj,me,rez,params,ndims,meta,sav)

rng(pi)


% m2q
for sessix = 1:numel(me)

    dt = params(sessix).dt;

    % transition data
    mdat = dat.mdat{sessix};
    nTrials = numel(mdat);


    % sort null/potent spaces by ve
    [~,nullix] = sort(rez(sessix).ve.null,'descend');
    nullix = nullix(1:ndims);

    [~,potentix] = sort(rez(sessix).ve.potent,'descend');
    potentix = potentix(1:ndims);

    nbouts = 5; % approximately max number of bouts in a trial
    nT = numel(obj(sessix).time);
    allnull = nan(nT,nTrials*nbouts,ndims);
    allpotent = nan(nT,nTrials*nbouts,ndims);
    allme = nan(nT,nTrials*nbouts,ndims);

    for dimix = 1:ndims

        for trix = 1:nTrials
            mdat_trix = mdat{trix};
            if isempty(mdat_trix)
                continue
            end

            sm = 11;
            null = mySmooth(rez(sessix).N_null(:,trix,nullix(dimix)),sm);
            potent = mySmooth(rez(sessix).N_potent(:,trix,potentix(dimix)),sm);
            move = mySmooth(me(sessix).move(:,trix),sm);
            medat = me(sessix).data(:,trix);

            % for each transition
            for i = 1:size(mdat_trix,2)
                tix = mdat_trix(1,i):mdat_trix(end,i);
                ts = tix .* dt;
                ts = ts - mean(ts);

                nullplot = null(tix);
                potentplot = potent(tix);
                meplot = medat(tix);

                % get center of tix
                dCenter = ceil(numel(tix)/2);
                center = round(nT / 2);
                try
                    allnull(center-dCenter:center+dCenter-1,(trix+i-1),dimix) = nullplot;
                    allpotent(center-dCenter:center+dCenter-1,(trix+i)-1,dimix) = potentplot;
                    allme(center-dCenter:center+dCenter-1,(trix+i)-1,dimix) = meplot;
                catch
                    allnull(center-dCenter:center+dCenter-2,(trix+i-1),dimix) = nullplot;
                    allpotent(center-dCenter:center+dCenter-2,(trix+i)-1,dimix) = potentplot;
                    allme(center-dCenter:center+dCenter-2,(trix+i)-1,dimix) = meplot;
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
    f.Position = [373    42   293   954];

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
    xlabel('Time to quiet (s)')
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
    xlabel('Time to quiet (s)')
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
    xlabel('Time to quiet (s)')
    ylabel('Null')
    ax(3).FontSize = 12;
    colorbar(ax(3)); %caxis([0 2.5])
    xline(0,'k:','LineWidth',2)
    colormap(cm)



    %%

    sgtitle([meta(sessix).anm ' - ' meta(sessix).date])

    if sav
        ppt.fig = f;
        ppt.figPath = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\figs\quiet-transitions';
        ppt.filename = 'quietToMoveSingleTrials';
        ppt.newVersion = 0;
        ppt.slideTitle = [meta(sessix).anm '  ' meta(sessix).date];
        myExportToPPTX(ppt)
    end


    clear ax



end



end