temprez = rez(end-3);
tempme = me(end-3);


trix = cell2mat(params(end-3).trialid(2:3)');

potent = temprez.N_potent;

potent = potent(:,trix,:);
potent = sum(potent.^2,3);


plotme = tempme.data(:,trix);

null = temprez.N_null;
null = null(:,trix,:);
null = sum(null.^2,3);

reward = obj(end-3).bp.ev.reward(trix);

%%

close all

f = figure;
f.Position = [680          66        1124         912];
dy = 75;

ax1 = subplot(1,3,1); hold on;
for i = 1:size(potent,2)
    plot(obj(1).time - (reward(i)-2.5) , potent(:,i) + i*dy, 'b')
end

ax = subplot(1,3,2); hold on;
for i = 1:size(potent,2)
    plot(obj(1).time , potent(:,i) + i*dy, 'r')
end

dy = 20;
ax = subplot(1,3,3); hold on;
for i = 1:size(plotme,2)
%     plot(obj(1).time  , plotme(:,i) + i*dy, 'k')
    plot(obj(1).time - (reward(i)-2.5) , plotme(:,i) + i*dy, 'k')
end



% align potent to reward, heatmap, gotta pad a matrix iwth nans or
% something



h = findobj(ax1,'Type','line');
for i = 1:numel(h)
    xs(:,i) = h(i).XData;
end

minmax = [min(min(xs)) max(max(xs))];

t = minmax(1):params(1).dt:minmax(2);

X = nan(numel(t), size(potent,2));

% now just put potent into X, aligned to it's specific xs












