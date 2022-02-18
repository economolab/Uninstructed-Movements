function lt = findLickTimesByCond(conds,params,lickStart,lickNum, plt)

[max_size, ~] = max(cellfun('size', params.trialid(conds), 1));
lt = nan(max_size,numel(conds));
for condix = 1:numel(conds) % for each condition
    c = conds(condix);
    trials = params.trialid{c};
    for trix = 1:numel(trials)
        trid = trials(trix);
        try
            lt(trix,condix) = lickStart{trid}(lickNum); % lickNum for current trial
        catch ME % catch when lickNum licks don't exist in current trial (expected)
            if (strcmp(ME.identifier,'MATLAB:badsubscript'))
                continue
            else
                rethrow(ME)
            end
        end
    end
end

% plot
fig = figure; hold on
for condix = 1:numel(conds)
        med = nanmedian(lt(:,condix));
        histogram(lt(:,condix),round(max_size/8), 'EdgeColor', 'none','FaceColor',plt.colors{condix})
        xline(med,'--','Color',plt.colors{condix},'HandleVisibility','off')
        plt.legend{condix} = [plt.legend{condix} ' ' num2str(med) ' (s)' ];
end
hold off
legend(plt.legend,'Location','best')
xlabel(['Lick ' num2str(lickNum) ' Start Time (s)'])
ylabel('# Trials')
title(['Lick ' num2str(lickNum)])
xlim([2 5])
ax = gca;
ax.FontSize = 20;

if plt.save
    savefn = ['Lick' num2str(lickNum) '_StartTimeDist'];
    saveas(fig,fullfile(pwd,savefn),'png')
end

end














