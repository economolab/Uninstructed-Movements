function plot_nonGC_moveTransitions_singleTrials_v3(dat,obj,me,rez,params,dim,meta)

rng(pi)

% m2q
for sessix = 1:numel(me)

    f = figure;
    f.Position = [20         595        1840         329];
    hold on;

    dt = params(sessix).dt;

    % transition data
    qdat = dat.qdat{sessix};

    % sort null/potent spaces by ve
    [~,nullix] = sort(rez(sessix).ve.null,'descend');
    nullix = nullix(dim);


    [~,potentix] = sort(rez(sessix).ve.potent,'descend');
    potentix = potentix(dim);

    nTrials = numel(qdat);

    nullplot = [];
    potentplot = [];
    meplot = [];
    xlineix = 0;

    for trix = 1:nTrials
        qdat_trix = qdat{trix};
        if isempty(qdat_trix)
            continue
        end

        npsm = 71;
        msm = 41;
        null = mySmooth(rez(sessix).N_null(:,trix,nullix),npsm);
        potent = mySmooth(rez(sessix).N_potent(:,trix,potentix),npsm);
        move = me(sessix).move(:,trix);
        medat = mySmooth(me(sessix).data(:,trix),msm);

        % for each transition
        for i = 1:size(qdat_trix,2)
            if rand > 0.15 % plot half of the transitions
                continue
            end

            tix = qdat_trix(1,i):qdat_trix(end,i);
            ix = find(tix == 0);
            if ~isempty(ix)
                tix(ix) = [];
            end
            ts = tix .* dt;
            ts = ts - mean(ts);

            nullplot = [nullplot ; null(tix)];

            potentplot = [potentplot ; potent(tix)];
            meplot = [meplot ; medat(tix)];
            xlineix = [xlineix ; numel(tix)];
            %             nrange(1) = min(min(nullplot),min(potentplot));
            %             nrange(2) = max(max(nullplot),max(potentplot));

            %             patchline(ts,normalize(meplot,'range',nrange),'EdgeColor','k','EdgeAlpha',0.35)
            %             patchline(ts,nullplot,'EdgeColor','r','EdgeAlpha',0.35)
            %             patchline(ts,potentplot,'EdgeColor','g','EdgeAlpha',0.35)



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
    ts = (1:numel(nullplot)) .* dt;


    nullplot = (nullplot .*2); % make it negative going
    potentplot = potentplot .*2;


    plot(ts,normalize(nullplot,'range',[0 1]),'Color',[62, 168, 105]./255,'LineWidth',2)
    plot(ts,normalize(potentplot,'range',[0 1])+0.9,'Color',[255, 56, 140]./255,'LineWidth',2)
    plot(ts,normalize(meplot,'range',[0 1])+0.65,'Color',[225, 144, 15]./255,'LineWidth',2)

    xlines = cumsum(xlineix) .* dt;
    xlines(1) = [];
    for i = 1:numel(xlines)
        xline(xlines(i),'k:','LineWidth',1)
    end
    xlim([0 ts(end)+0.5])
    title([meta(sessix).anm ' - ' meta(sessix).date])
    xlabel('seconds (quiet to move transitions)')



end



end