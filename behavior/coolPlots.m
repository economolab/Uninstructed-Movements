clear,clc,close all


datapth = '/Users/Munib/Documents/Economo-Lab/code/data';
datafn = 'data_structure_JGR2_2021-11-16.mat';

load(fullfile(datapth,datafn))

% load('Y:\eng_research_economo\JEB\Experiments\JGR2\Analysis\2021-11-16\data_structure_JGR2_2021-11-16.mat');

traj = obj.traj{1}; 

% % features
% % 1 - tongue
% % 2 - left_tongue
% % 3 - right_tongue
% % 4 - jaw
% % 5 - trident
% % 6 - nose
% % 7 - lickport

%%


% % get midpoint bw jaw and trident across all trials
% mx = [];
% my = [];
% for trix = 1:numel(traj)
%     jx = traj(trix).ts(:,1,4);
%     jy = traj(trix).ts(:,2,4);
%     tx = traj(trix).ts(:,1,5);
%     ty = traj(trix).ts(:,2,5);
%     midx = (jx + tx)./2;
%     midy = (jy + ty)./2;
% 
%     mx = [mx ; midx];
%     my = [my ; midy];
% end
% medx = nanmedian(mx); % 167.64
% medy = nanmedian(my); % 142.57

medx = 167.64;
medy = 142.57;

% concatenate displacement from (medx,medy) 
% l = obj.bp.L&obj.bp.hit;
% r = obj.bp.R&obj.bp.hit;
l = obj.bp.L&obj.bp.hit&~obj.bp.autowater&~obj.bp.early;
r = obj.bp.R&obj.bp.hit&~obj.bp.autowater&~obj.bp.early;
le = obj.bp.L&obj.bp.miss&~obj.bp.autowater&~obj.bp.early;
re = obj.bp.R&obj.bp.miss&~obj.bp.autowater&~obj.bp.early;
law = obj.bp.L&obj.bp.hit&obj.bp.autowater&~obj.bp.early;
raw = obj.bp.R&obj.bp.hit&obj.bp.autowater&~obj.bp.early;

dl = []; % left
dr = []; % right
dre = []; % right error
dle = []; % left error
draw = []; % right autowater
dlaw = []; % left autowater
for trix = 1:obj.bp.Ntrials%numel(traj)

    tm = (1:size(traj(trix).ts,1))./400;
    jx = traj(trix).ts(:,1,4);
    jy = traj(trix).ts(:,2,4);
    tx = traj(trix).ts(:,1,5);
    ty = traj(trix).ts(:,2,5);
    midx = (jx + tx)./2;
    midy = (jy + ty)./2;    

    temp = sqrt((midx-medx).^2 + (midy-medy).^2); % displacement from median (medianx and median y hardcoded)

    if l(trix)
        dl = [dl ; temp];
    elseif r(trix)
        dr = [dr ; temp];
    elseif le(trix)
        dle = [dle ; temp];
    elseif re(trix)
        dre = [dre ; temp];
    elseif law(trix)
        dlaw = [dlaw; temp];
    elseif raw(trix)
        draw = [draw; temp];
    end
end

figure; 
subplot(411)
plot(dl,'r'); title('Left 2afc')
subplot(412)
plot(dr,'b'); title('Right 2afc')
subplot(413)
plot(dlaw,'Color',[1 0.5 0.5]); title('Left aw')
subplot(414); 
plot(draw,'Color',[0.5 0.5 1]); title('Right aw')


%%

figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2); 

for trix = 150:numel(traj)
    tm = (1:size(traj(trix).ts,1))./400;
    jx = traj(trix).ts(:,1,4);
    jy = traj(trix).ts(:,2,4);
    tx = traj(trix).ts(:,1,5);
    ty = traj(trix).ts(:,2,5);
    midx = (jx + tx)./2;
    midy = (jy + ty)./2;

    d = sqrt((midx-166).^2 + (midy-142).^2); % displacement from median (medianx and median y hardcoded)
    plot(ax1, tm,d,'k-'); xlabel('time'); ylabel('mid jaw trident displacement');
    scatter3(ax2, midx, -midy, 10e-2.*medfilt1(myDiff(d,1/400), 5), 50, tm, '.'); colormap(jet); %axis image
    xlabel('x'); ylabel('y'); zlabel('velocity');
    colorbar;
    pause
    cla(ax1);
    cla(ax2);
end

%% do above for for random subset of 20 aw and 20 2afc trials
k = 20;

afc = obj.bp.hit&~obj.bp.autowater&~obj.bp.early;
aw = obj.bp.hit&obj.bp.autowater&~obj.bp.early;

trials.afc = randsample(find(afc),k);
trials.aw = randsample(find(aw),k);

f(1) = figure(1);
ax(1) = axes(f(1));
set(f(1),'Position',[-1809         193         741         517]);
f(2) = figure(2);
ax(2) = axes(f(2));
set(f(2),'Position',[-1061         195         842         513]);


hold(ax(1),'on')
hold(ax(2),'on')

for trix = 1:k
    curtrix = trials.afc(trix);
    tm = (1:size(traj(curtrix).ts,1))./400;
%     ix = find(tm>obj.bp.ev.goCue(curtrix),1,'first');
    ix = numel(tm);
    tm=tm(1:ix);
    jx = traj(curtrix).ts(1:ix,1,4);
    jy = traj(curtrix).ts(1:ix,2,4);
    tx = traj(curtrix).ts(1:ix,1,5);
    ty = traj(curtrix).ts(1:ix,2,5);
    midx = (jx + tx)./2;
    midy = (jy + ty)./2;
    d = sqrt((midx-166).^2 + (midy-142).^2); % displacement from median (medianx and median y hardcoded)
    scatter3(ax(1), midx, -midy, 10e-2.*medfilt1(myDiff(d,1/400), 5), 50, tm, '.'); 
    colormap(ax(1),jet); %axis image
    xlabel(ax(1),'x'); ylabel(ax(1),'y'); zlabel(ax(1),'velocity');
    colorbar(ax(1));
    ax(1).FontSize = 30;
    ax(1).XTick = [];
    ax(1).YTick = [];
    ax(1).ZTick = [];
    
    curtrix = trials.aw(trix);
    tm = (1:size(traj(curtrix).ts,1))./400;
    ix = find(tm>obj.bp.ev.goCue(curtrix),1,'first');
    tm=tm(1:ix);
    jx = traj(curtrix).ts(1:ix,1,4);
    jy = traj(curtrix).ts(1:ix,2,4);
    tx = traj(curtrix).ts(1:ix,1,5);
    ty = traj(curtrix).ts(1:ix,2,5);
    midx = (jx + tx)./2;
    midy = (jy + ty)./2;
    d = sqrt((midx-166).^2 + (midy-142).^2); % displacement from median (medianx and median y hardcoded)
    scatter3(ax(2), midx, -midy, 10e-2.*medfilt1(myDiff(d,1/400), 5), 50, tm, '.'); 
    colormap(ax(2),jet); %axis image
    xlabel(ax(2),'x'); ylabel(ax(2),'y'); zlabel(ax(2),'velocity');
    colorbar(ax(2));
    
end

%%
featNames = traj(1).featNames;

legendString = {'tongueX', 'tongueY', 'leftTongueX','leftTongueY',...
    'rightTongueX','rightTongueY','jawX','jawY','tridentX','tridentY','noseX','noseY',...
    'lickPortX','lickPortY'};




clr = {[0, 0.4470, 0.7410],
    [0.8500, 0.3250, 0.0980],
    [0.9290, 0.6940, 0.1250],
    [0.4940, 0.1840, 0.5560],
    [0.4660, 0.6740, 0.1880],
    [0.3010, 0.7450, 0.9330],
    [0.6350, 0.0780, 0.1840]};

ls = {'--','-'};


close all
fig = figure; 
fig.Position = [1          41        2560        1323];
for trix = 1:numel(traj)
    hold on;
    for fix = 1:numel(featNames)
        tm = (1:size(traj(trix).ts,1))./400;
        for cix = 1:2 % coord
            ts = traj(trix).ts(:,cix,fix);
            plot(tm,ts,ls{cix},'Color',clr{fix})
        end
    end
    title(['Trial ' num2str(trix) ' ' traj(trix).fn(1:51)],'Interpreter','none')
    xlabel('Time (s)')
    ylabel('x (dashed)  ||   y (solid)')
    legend(legendString,'location','bestoutside')
    ax = fig.CurrentAxes;
    ax.FontSize = 20;

    pause
    clf
end