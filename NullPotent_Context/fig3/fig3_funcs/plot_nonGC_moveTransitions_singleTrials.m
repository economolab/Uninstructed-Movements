function plot_nonGC_moveTransitions_singleTrials(dat,obj,me,rez,params,dim,meta)

rng(pi)

% m2q
for sessix = 1:numel(me)
    f = figure;
    f.Position = [403         163        1106         792];
    hold on;

    dt = params(sessix).dt;

    % transition data
    mdat = dat.mdat{sessix};

    % sort null/potent spaces by ve
    [~,nullix] = sort(rez(sessix).ve.null,'descend');
    nullix = nullix(dim);


    [~,potentix] = sort(rez(sessix).ve.potent,'descend');
    potentix = potentix(dim);

    nTrials = numel(mdat);

    for trix = 1:nTrials
        mdat_trix = mdat{trix};
        if isempty(mdat_trix)
            continue
        end

        sm = 41;
        null = mySmooth(rez(sessix).N_null(:,trix,nullix),sm);
        potent = mySmooth(rez(sessix).N_potent(:,trix,potentix),sm);
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

            patchline(ts,normalize(meplot,'range',[-1 1]) + 3,'EdgeColor','k','EdgeAlpha',0.35)
            patchline(ts,nullplot,'EdgeColor','r','EdgeAlpha',0.35)
            patchline(ts,potentplot+5,'EdgeColor','g','EdgeAlpha',0.35)



            %             plot(ts,normalize(medat(tix),'range',[-1 1]) + 5,'k')
            %             plot(ts,null(tix),'r')
            %             plot(ts,potent(tix)+10,'g')

            %             hold on;
            %             plot(ts,normalize(medat(tix),'range',[1.5 2.5]),'k')
            %             plot(ts,normalize(null(tix),'range',[0 1]),'r')
            %             plot(ts,normalize(potent(tix),'range',[3 4]),'g')
            %             pause
            %             hold off
            %             clf
        end

    end
    title([meta(sessix).anm ' - ' meta(sessix).date])
    xlabel('Time to quiet (s)')

    %     % mean
    %     meannull = nanmean(allnull,2);
    %     meanpotent = nanmean(allpotent,2);
    %     meanme = normalize(nanmean(allme,2),'range',[0 1]);
    %     f = figure;
    %     hold on
    %     plot(obj(sessix).time,meanme,'k','LineWidth',2)
    %     plot(obj(sessix).time,meannull,'r','LineWidth',2)
    %     plot(obj(sessix).time,meanpotent,'g','LineWidth',2)
    %     xlim([-1 1])
    %     title([meta(sessix).anm ' - ' meta(sessix).date])
    %     xlabel('Time to quiet (s)')


end



end