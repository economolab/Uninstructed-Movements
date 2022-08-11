function varexpPlots(rez)

null_total = zeros(size(rez));
potent_total = zeros(size(rez));
null_prep = zeros(size(rez));
null_move = zeros(size(rez));
potent_move = zeros(size(rez));
potent_prep = zeros(size(rez));
for i = 1:numel(rez)
    null_total(i) = rez(i).ve.null_total;
    potent_total(i) = rez(i).ve.potent_total;
    null_prep(i) = rez(i).ve.null_prep;
    null_move(i) = rez(i).ve.null_move;
    potent_move(i) = rez(i).ve.potent_move;
    potent_prep(i) = rez(i).ve.potent_prep;
end

violincols = [50, 168, 82; 168, 50, 142] ./ 225;
varexp_full = [null_total ; potent_total]';
f = figure; ax = axes(f);
vs = violinplot(varexp_full,{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols);
ylabel('Normalized Variance Explained (Whole Trial)')
ylim([0,1])
ax = gca;
ax.FontSize = 20;
% 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace';
% fn = 've_total';
% mysavefig(f,pth,fn);

violincols = [50, 168, 82; 168, 50, 142; 168, 50, 142; 50, 168, 82] ./ 225;
varexp_sep = [null_prep; null_move ; potent_move; potent_prep]';
varexp_sep(:,1) = varexp_sep(:,1) - 0.001*rand(numel(rez),1);
f = figure; ax = axes(f);
vs = violinplot(varexp_sep,{'Null, Non-move','Null, Move', 'Potent,Move','Potent, Non-move'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1},  'ViolinColor', violincols);
ylabel('Normalized Variance Explained')
ylim([0,1])
ax = gca;
ax.FontSize = 20;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/pcaNullSpace';
% fn = 've_allsep';
% mysavefig(f,pth,fn);

end