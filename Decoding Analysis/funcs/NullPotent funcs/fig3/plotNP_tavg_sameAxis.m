close all

cols = getColors;
% fns = fieldnames(cols);
% for i = 1:numel(fns)
%     clrs{i} = cols.(fns{i});
% end
clrs{2} = [99, 3, 3] ./ 255;
clrs{1} = [245, 91, 112] ./ 255; %'k'; %
clrs{4} = [3, 6, 99] ./ 255;
clrs{3} = [101, 88, 245] ./ 255;

% clrs{1} = cols.rhit_aw;
% clrs{2} = cols.rhit;
% clrs{3} = cols.lhit_aw;
% clrs{4} = cols.lhit;

for sessix = 1:numel(meta)
    trix{1} = params(sessix).trialid{2};
    trix{2} = params(sessix).trialid{3};

    %     f = figure;
    %     ax = subplot(2,1,1);
    %     hold on
    %     for i = 1:numel(trix)
    %         null{i} = rez(sessix).N_null(:,trix{i},:);
    %         null{i} = sum(null{i}.^2,3);
    %         nullmean = mean(null{i},2);
    %         nullstd = std(null{i},[],2) ./ sqrt(size(null{i},2));
    %         shadedErrorBar(obj(sessix).time,nullmean,nullstd,{'Color',clrs{i},'LineWidth',2},0.2,ax)
    %     end
    %
    %     ax = subplot(2,1,2);
    %     hold on
    %     for i = 1:numel(trix)
    %         potent{i} = rez(sessix).N_potent(:,trix{i},:);
    %         potent{i} = sum(potent{i}.^2,3);
    %         potentmean = mean(potent{i},2);
    %         potentstd = std(potent{i},[],2) ./ sqrt(size(potent{i},2));
    %         shadedErrorBar(obj(sessix).time,potentmean,potentstd,{'Color',clrs{i},'LineWidth',2},0.2,ax)
    %     end

    f = figure;
    ax = gca;
    hold on
    cols2use = {[1 2],[3 4]};
    for i = 1:numel(trix)
        null{i} = rez(sessix).N_null(:,trix{i},:);
        null{i} = sum(null{i}.^2,3);
        nullmean = mean(null{i},2);
        nullstd = std(null{i},[],2) ./ sqrt(size(null{i},2));
        shadedErrorBar(obj(sessix).time,nullmean,nullstd,{'Color',clrs{cols2use{i}(1)},'LineWidth',2},0.2,ax)

        potent{i} = rez(sessix).N_potent(:,trix{i},:);
        potent{i} = sum(potent{i}.^2,3);
        potentmean = mean(potent{i},2);
        potentstd = std(potent{i},[],2) ./ sqrt(size(potent{i},2));
        shadedErrorBar(obj(sessix).time,potentmean,potentstd,{'Color',clrs{cols2use{i}(2)},'LineWidth',2},0.2,ax)
    end

    sample = mode(obj(sessix).bp.ev.sample) - 2.5;
    delay = mode(obj(sessix).bp.ev.delay) - 2.5;
    xline(0,'k--')
    xline(sample,'k:')
    xline(delay,'k:')
    xlim([obj(sessix).time(5) 2])
%     ylim([0 250])

end