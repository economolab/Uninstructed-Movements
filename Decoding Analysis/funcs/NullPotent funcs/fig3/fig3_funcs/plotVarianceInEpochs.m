function plotVarianceInEpochs(rez,me,params,cond2use)
%%
maxdims = 10; % only get variance in top 10 variance explaining dims

types = {'null-move','null-prep','potent-move','potent-prep'};
vars = nan(numel(rez),4,maxdims); % (sessions,ntypes,maxdims)

for sessix = 1:numel(rez)

    % trials
    trix = cell2mat(params(sessix).trialid(cond2use)');
    
    % move time points
    move = me(sessix).move(:,trix);

    % -null
    % sort dimensions by variance explained (most to least)
    [~,ix] = sort(rez(sessix).ve.null,'descend');
    temp = rez(sessix).N_null(:,trix,ix);
    temp2 = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
    for dimix = 1:min(maxdims,size(temp,3))
        temp3 = temp2(move(:),dimix);
        vars(sessix,1,dimix) = var(temp3);
        temp3 = temp2(~move(:),dimix);
        vars(sessix,2,dimix) = var(temp3);
    end
    % -potent
    % sort dimensions by variance explained (most to least)
    [~,ix] = sort(rez(sessix).ve.potent,'descend');
    temp = rez(sessix).N_potent(:,trix,ix);
    temp2 = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
    for dimix = 1:min(maxdims,size(temp,3))
        temp3 = temp2(move(:),dimix);
        vars(sessix,3,dimix) = var(temp3);
        temp3 = temp2(~move(:),dimix);
        vars(sessix,4,dimix) = var(temp3);
    end
end


%% plot

cols = linspecer(10);

f=figure; hold on;
ax = gca;
div = 1;
xs = [1 2 5 6];
for i = 1:size(vars,2)
    barplot = mean(vars(:,i,:),3); % mean across dims
    h(i) = bar(xs(i),mean(barplot)); % mean across dims and sessions
    h(i).FaceColor = 'k';
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    for dimix = 1:maxdims
        scatterplot = mean(vars(:,i,dimix)); % mean variance for a type and dim
        scatter(xs(i)+(0.2*rand),scatterplot,60,'MarkerFaceColor',cols(dimix,:)./div, ...
            'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',1, ...
            'MarkerFaceAlpha',0.7)
    end
%     errorbar(h(i).XEndPoints,mean(vars(:,i)),std(vars(:,i))./sqrt(size(vars,1)),'LineStyle','none','Color','k','LineWidth',1);
end
ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xticklabels(types);

ylabel('Mean variance across trials and sessions')
ax.FontSize = 12;







%%%%%%%%%%%%%%%%

% var_.null_move = cell(numel(rez),1);
% var_.null_prep = cell(numel(rez),1);
% var_.potent_move = cell(numel(rez),1);
% var_.potent_prep = cell(numel(rez),1);
% for sessix = 1:numel(rez)
%     % trials
%     trix = [];
%     for condix = 1:numel(cond2use)
%         trix = [trix ; params(sessix).trialid{cond2use(condix)}];
%     end
% 
%     move = me(sessix).move(:,trix);
% 
% 
%     % -null
%     % sort dimensions by variance explained (most to least)
%     [~,ix] = sort(rez(sessix).ve.null,'descend');
%     temp = rez(sessix).N_null(:,trix,ix);
%     temp2 = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
%     for dimix = 1:size(temp,3)
%         temp3 = temp2(move(:),dimix);
%         var_.null_move{sessix}(dimix) = var(temp3);
%         temp3 = temp2(~move(:),dimix);
%         var_.null_prep{sessix}(dimix) = var(temp3);
%     end
%     % -potent
%     % sort dimensions by variance explained (most to least)
%     [~,ix] = sort(rez(sessix).ve.potent,'descend');
%     temp = rez(sessix).N_potent(:,trix,ix);
%     temp2 = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
%     for dimix = 1:size(temp,3)
%         temp3 = temp2(move(:),dimix);
%         var_.potent_move{sessix}(dimix) = var(temp3);
%         temp3 = temp2(~move(:),dimix);
%         var_.potent_prep{sessix}(dimix) = var(temp3);
%     end
% end
% 
% vars(:,1) = mean(cell2mat(var_.null_prep),2); % null_prep
% vars(:,2) = mean(cell2mat(var_.null_move),2); % null_move
% vars(:,3) = mean(cell2mat(var_.potent_move),2); % potent_move
% vars(:,4) = mean(cell2mat(var_.potent_prep),2); % null_move
% 
% cols(1,:) = [50, 191, 83];
% cols(2,:) = [255, 128, 48];
% cols(3,:) = [232, 53, 226];
% cols(4,:) = [53, 226, 232];
% cols = cols./255;

% f=figure; hold on;
% ax = gca;
% div = 1;
% xs = [1 2 5 6];
% for i = 1:size(vars,2)
%     h(i) = bar(xs(i),mean(vars(:,i)));
%     h(i).FaceColor = cols(i,:);
%     h(i).EdgeColor = 'none';
%     h(i).FaceAlpha = 0.5;
%     scatter(xs(i)*ones(size(vars(:,i))),vars(:,i),60,'MarkerFaceColor',cols(i,:)./div, ...
%             'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.35, ...
%             'MarkerFaceAlpha',0.7)
%     errorbar(h(i).XEndPoints,mean(vars(:,i)),std(vars(:,i))./sqrt(size(vars,1)),'LineStyle','none','Color','k','LineWidth',1);
% end
% 
% ylim([-0.001 ax.YLim(2)])
% ax.XTick = xs;
% xticklabels({'null prep','null move','potent move','potent prep'})
% ylabel('Mean variance across trials and dims')
% ax.FontSize = 14;


end