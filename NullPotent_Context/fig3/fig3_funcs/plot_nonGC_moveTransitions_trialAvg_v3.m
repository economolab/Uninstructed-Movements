function plot_nonGC_moveTransitions_trialAvg_v3(dat,obj,me,rez,params,meta,ndims)

rng(pi)

% generate colors
h = 0.0;
s = linspace(0.4,1,ndims);
v = linspace(0.5,1,ndims);
for i = 1:ndims
    hsv(i,:) = [h s(i) v(i)];
end
rgb.null = hsv2rgb(hsv);

h = 0.38;
for i = 1:ndims
    hsv(i,:) = [h s(i) v(i)];
end
rgb.potent = hsv2rgb(hsv);


% m2q
for sessix = 1:numel(me)
    f = figure;
    f.Position = [403         163        1106         792];
    ax = gca;
    hold on

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
    allnull = nan(500,nTrials*nbouts,ndims);
    allpotent = nan(500,nTrials*nbouts,ndims);
    allme = nan(500,nTrials*nbouts,ndims);

    for dimix = 1:ndims

        for trix = 1:nTrials
            mdat_trix = mdat{trix};
            if isempty(mdat_trix)
                continue
            end

            sm = 41;
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
    % mean
    meannull = squeeze(nanmean(allnull,2));
    meanpotent = squeeze(nanmean(allpotent,2));
    meanme = squeeze(normalize(nanmean(allme,2),'range',[-0.3 0.3]));
    meanme = meanme(:,1);

    nullmeanplot = mean(meannull,2);
    nullerrplot = std(meannull,[],2) ./ sqrt(ndims);
    potentmeanplot = mean(meanpotent,2);
    potenterrplot = std(meanpotent,[],2) ./ sqrt(ndims);


    plot(obj(sessix).time,meanme,'m','LineWidth',3)
    shadedErrorBar(obj(sessix).time,nullmeanplot,nullerrplot,{'Color','r','LineWidth',2},0.5,ax)
    shadedErrorBar(obj(sessix).time,potentmeanplot,potenterrplot,{'Color','g','LineWidth',2},0.5,ax)

    xline(0,'k:','LineWidth',2)

    xlim([-1 1])
    title([meta(sessix).anm ' - ' meta(sessix).date])
    xlabel('Time to quiet (s)')

end



end