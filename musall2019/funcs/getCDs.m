function getCDs(obj,dat,params,fit_type)


ev.sample = obj.bp.ev.sample;
ev.delay = obj.bp.ev.delay;
ev.goCue = obj.bp.ev.goCue;
ev.(params.alignEvent) = obj.bp.ev.(params.alignEvent);

rez.time = obj.time;
rez.psth = dat;
% rez.psth = standardizePSTH(obj);
rez.condition = params.condition;
rez.alignEvent = params.alignEvent;
rez.ev = ev;


% cd early mode

e1 = mode(ev.delay) - 0.5 - mode(ev.(params.alignEvent));
e2 = mode(ev.delay) - 0.1 - mode(ev.(params.alignEvent));
%     e1 = mode(ev.sample) + 0.4 - mode(ev.(params.alignEvent));
%     e2 = mode(ev.sample) + 0.8 - mode(ev.(params.alignEvent));

times.early = rez.time>e1 & rez.time<e2;
rez.cdEarly_mode = calcCD(rez,times.early);


% cd late mode

e1 = mode(ev.goCue) - 0.5 - mode(ev.(params.alignEvent));
e2 = mode(ev.goCue) - 0.1 - mode(ev.(params.alignEvent));

times.late = rez.time>e1 & rez.time<e2;
rez.cdLate_mode = calcCD(rez,times.late);


% cd go mode

e1 = mode(ev.goCue) + 0.02 - mode(ev.(params.alignEvent));
e2 = mode(ev.goCue) + 0.42 - mode(ev.(params.alignEvent));

times.go = rez.time>e1 & rez.time<e2;
rez.cdGo_mode = calcCD(rez,times.go);

% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(rez.cdLate_mode),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end


% projections and normalize (not normalizing currently since I
% standardize PSTHs beforehand)
% when pooling trajectories across sessions from hidehikos ppn paper:
% CD_late projections normalized by mean activity just before go cue
% (-0.1<t<t_go)
% CD_go projections normalized by mean activity after go cue
% (t_go<t<0.4)

%     normTimes{1} = rez.time>e1 & rez.time<e2; % sample
%     normTimes{2} = rez.time>-0.4 & rez.time<0; % delay
%     normTimes{3} = rez.time>0 & rez.time<0.4; % go

cond = [1 2];
for i = 1:numel(fns)
    tempmode = rez.(fns{i});
    for j = 1:numel(cond)
        c = cond(j);

        tempdat = rez.psth(:,:,c)*rez.(fns{i});

        %             normfactor = abs(nanmean(tempdat(normTimes{i})));
        normfactor = 1;
        %             normfactor = max(tempdat);

        rez.([fns{i}(1:end-5) '_latent'])(:,j) = tempdat ./ normfactor;
    end
end

clear cond

% variance explained

for i = 1:numel(fns)
    psth = rez.psth;
    datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
    datacov(isnan(datacov)) = 0;
    eigsum = sum(eig(datacov));
    rez.varexp.(fns{i}(1:end-5)) = var_proj(rez.(fns{i}), datacov, eigsum);
end


% concatenate latents from all sessions into a (time,sessions) matrix, find mean and stderror across sessions

fns = patternMatchCellArray(fieldnames(rez(1)),{'latent'},'all');

for fnix = 1:numel(fns)
    fn = fns{fnix};
    % preallocate
    cd.(fn){1} = [];
    cd.(fn){2} = [];
    % concatenate
    for i = 1:numel(rez)
        for j = 1:2
            cd.(fn){j} = [cd.(fn){j} rez(i).(fn)(:,j)];
        end
    end
    % mean and stderror across session (just for plotting)[p
    cd.mean.(fn) = [nanmean(cd.(fn){1},2) nanmean(cd.(fn){2},2)];
    cd.stderr.(fn) = [nanstd(cd.(fn){1},[],2) nanstd(cd.(fn){2},[],2)] ./ numel(rez); % std error
end

% average varexp from all sessions
fns = fieldnames(rez(1).varexp);
for i = 1:numel(rez)
    for j = 1:numel(fns)
        vexp(i,j) = rez(i).varexp.(fns{j});
    end
end
mean_ve = mean(vexp);



% selectivity
for i = 1:numel(rez)
    for j = 1:numel(fns)
        selectivity{i,j} = cd.mean.([fns{j} '_latent'])(:,1) - cd.mean.([fns{j} '_latent'])(:,2);
    end
end

%% plot projections onto coding directions
clrs = getColors();
lw = 3;
alph = 0.5;

sample = mode(rez(1).ev.sample - rez(1).ev.(params.alignEvent));
delay = mode(rez(1).ev.delay - rez(1).ev.(params.alignEvent));

fns = fieldnames(cd.mean);

sav = 0;
f = figure;
f.Position = [318   679   995   281];
t = tiledlayout('flow');
for i = 1:numel(fns)
    %     f = figure; hold on
    %     ax = gca;
    ax = nexttile; hold on;
    tempmean = cd.mean.(fns{i});
    temperror = cd.stderr.(fns{i});
    shadedErrorBar(rez(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(rez(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)

    xlim([rez(1).time(10);rez(1).time(end)])
    %     ylims = [min(min(tempmean))-5, max(max(tempmean))+5];
    %     ylim(ylims);

%     title(['$ CD_{' lower(fns{i}(3:end-7)) '}$'  ' | MeanVE=' num2str(round(mean_ve(i),2)) ],'Interpreter','latex')
    title(['$ CD_{' lower(fns{i}(3:end-7)) '}$'],'Interpreter','latex')
    %     xlabel('Time (s) from go cue')
    %     ylabel('Activity (a.u.)')
    %     ax = gca;
    ax.FontSize = 15;
%     ax.XTick = [];
    %     ax.YTick = [];

    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)

    curmodename = fns{i};
    timefieldname = [lower(curmodename(3:end-7))];
    shadetimes = obj.time(times.(timefieldname));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
        y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%     y = [-60 -60 50 50];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';

%         ylim(ylims);
%     ylim([-60 50]);


    %     if sav
    %         pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/cd';
    %         fn = [fns{i} '_anmList1_sessionList1_w_excludedsessions_sm_' num2str(params.smooth)];
    %         mysavefig(f(i),pth,fn);
    %     end

    hold off

end
title(t,fit_type,'FontSize',15)
xlabel(t,'Time (s) from go cue','FontSize',15)
ylabel(t, 'au','FontSize',15)


%% selectivity

f = figure;
f.Position = [469   310   995   281];
t = tiledlayout('flow');
for i = 1:numel(selectivity)
  %     f = figure; hold on
    %     ax = gca;
    ax = nexttile; hold on;
    plot(rez(1).time,selectivity{i},'Color',[211, 89, 227]./255,'LineWidth',lw)
    zline = yline(0,'--','LineWidth',2);
    zline.Color = [227, 89, 89]./255;

    xlim([rez(1).time(10);rez(1).time(end)])

    title(['$ CD_{' lower(fns{i}(3:end-7)) '}$'],'Interpreter','latex')

    ax.FontSize = 15;

    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)

    curmodename = fns{i};
    timefieldname = [lower(curmodename(3:end-7))];
    shadetimes = obj.time(times.(timefieldname));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
        y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';



    hold off
end
title(t,fit_type,'FontSize',15)
xlabel(t,'Time (s) from go cue','FontSize',15)
ylabel(t, 'Selctivity (a.u.)','FontSize',15)

%% vexp
% violincols = [50, 168, 82; 168, 50, 142] ./ 225;
% f = figure; ax = axes(f);
% 
% tempvexp = [vexp sum(vexp,2)];
% vfns = fieldnames(rez(1).varexp);
% vfns{end+1} = 'sum';
% vs = violinplot(tempvexp,vfns,...
%     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1});% , 'ViolinColor', violincols);
% ylabel('Fraction of VE')
% ylim([0,1])
% ax = gca;
% ax.FontSize = 20;


%% orthogonality of modes
% close all
% dotProductModes(rez,orthModes,'Orthogonality of Coding Directions')
% f=gcf;
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/cd';
% fn = ['orthogonality_anmList1_sessionList1_w_excludedsessions_sm_' num2str(params.smooth)];
% mysavefig(f,pth,fn);


%% Helper Functions


function cd = calcCD(rez,times)
tempdat = rez.psth(:,:,[1,2]);
mu = squeeze(mean(tempdat(times,:,:),1));
%     toswitch = find(mu(:,1) < mu(:,2));
%     for i = 1:numel(toswitch)
%         new = flip(mu(toswitch(i),:));
%         mu(toswitch(i),:) = new;
%     end
sd = squeeze(std(tempdat(times,:,:),[],1));
%     for i = 1:numel(toswitch)
%         new = flip(sd(toswitch(i),:));
%         sd(toswitch(i),:) = new;
%     end
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)





