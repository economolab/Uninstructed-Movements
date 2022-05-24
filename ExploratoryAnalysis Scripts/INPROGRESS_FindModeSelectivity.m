%%
% Define time epoch that you want to assess selectivity for
start_ix = find(rez.time<-1.1975 & rez.time>-1.2075);     % 1.2025 s before first lick
stop_ix = find(rez.time<-0.2975 & rez.time>-0.3075);      % to 0.3025 s before first lick

rez.selectivity = NaN(1,numel(objs));
rez.normselectivity = NaN(1,numel(objs));
removeEarly.selectivity = NaN(1,numel(objs));
removeRandom.selectivity = NaN(1,numel(objs));
for i = 1:numel(objs)
    % Find choice mode selectivity in normal modes
    selectivity = rez.latentChoice{i,1} - rez.latentChoice{i,2};     % Take the difference in projections onto normal modes for R and L
    if ~isempty(selectivity)
        avg_selectivity = mean(selectivity(start_ix:stop_ix),'omitnan');        % Find the average selectivity during the specified time epoch
        rez.selectivity(i) = abs(avg_selectivity);               % Take absolute value of selectivity (only care about magnitude)
        m1 = max(abs(rez.latentChoice{i,1}));
        m2 = max(abs(rez.latentChoice{i,2}));
        m = max([m1 m2]);
        rez.normselectivity(i) = rez.selectivity(i)/m;
    end

    % Find choice mode selectivity in modes defined without early movements
    selectivity = removeEarly.latentChoice{i,1} - removeEarly.latentChoice{i,2};     % Take the difference in projections onto earlyRemove modes for R and L
    if ~isempty(selectivity)
        avg_selectivity = mean(selectivity(start_ix:stop_ix),"omitnan");
        removeEarly.selectivity(i) = abs(avg_selectivity);
    end

    % Find choice mode selectivity in modes defined with random trials
    % removed
    selectivity = removeRandom.latentChoice{i,1} - removeRandom.latentChoice{i,2};     % Take the difference in projections onto earlyRemove modes for R and L
    if ~isempty(selectivity)
        avg_selectivity = mean(selectivity(start_ix:stop_ix),'omitnan');
        removeRandom.selectivity(i) = abs(avg_selectivity);
    end
end
%%
% Find ratio of selectivity in choice mode when early movements are removed
% to the normal choice mode
% < 1 = selectivity was reduced by removing early move trials
% 1 = selectivity was unchanged 
% > 1 = selectivity was increased by removing early move trials
x = removeEarly.selectivity./rez.selectivity;  
r = removeRandom.selectivity./rez.selectivity;
p = NaN(1,numel(objs));
for ii = 1:length(x)
    obj = objs{ii};
    if isfield(obj,'pctEarlyMove')
        p(ii) = obj.pctEarlyMove;
    end
    if isnan(x(ii))
        r(ii) = NaN;
        p(ii) = NaN;
        rez.normselectivity(ii) = NaN;
    end
end

figure(1);
edges = linspace(min(rez.normselectivity),max(rez.normselectivity),20);
subplot(3,1,1)
histogram(rez.normselectivity,edges)
ylabel('# of sessions')
xlabel('R - L selectivity')
title('All trials included')

subplot(3,1,2)
edges = linspace(min(x),max(x),20);
histogram(x,edges)
ylabel('# of sessions')
xlabel('Fraction of selectivity remaining')
title('Early move trials removed')

subplot(3,1,3)  
edges = linspace(min(r),max(r),20);
histogram(r,edges)
ylabel('# of sessions')
xlabel('Fraction of selectivity remaining')
title('Random trials removed')

figure(2)
histogram(x,edges)
hold on
histogram(r,edges)
ylabel('# of sessions')
xlabel('Fraction of selectivity remaining')
legend('Early move trials','Random trials')


