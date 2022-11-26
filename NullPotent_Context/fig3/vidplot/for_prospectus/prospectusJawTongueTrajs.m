% can run this script after running st_elsayed/pca/kaufman.m

% set session data to use
anm = 'JEB15';
date = '2022-07-27';

anms = {meta(:).anm};
dates = {meta(:).date};

anmix = ismember(anms,anm);
dateix = ismember(dates,date);
sessix = find(all([anmix;dateix],1));

dat.rez = rez(sessix);
dat.obj = obj(sessix);
dat.params = params(sessix);
dat.meta = meta(sessix);
dat.me = me(sessix);

%%

close all

rng(pi+1)
trix = find(~dat.obj.bp.early&(dat.obj.bp.hit|dat.obj.bp.miss));
trix = randsample(trix,10,false);

gocue = 2.5;
sample = mode(dat.obj.bp.ev.sample) - gocue;
delay = mode(dat.obj.bp.ev.delay) - gocue;

cols = linspecer(3);
f = figure; hold on
dy = 1;
for itrix = 1:numel(trix)

    % get obj trajectories
    view = 1; % side cam
    trial = trix(itrix);

    frametimes = dat.obj.traj{view}(trial).frameTimes - 0.5;



    % resmaple obj.time to frametimes, but first need to get frametimes from
    % -2.5 from go cue to 2.5 to go cue
    tmin = dat.params.tmin;
    tmax = dat.params.tmax;

    [~,ix1] = min(abs(frametimes - gocue - tmin));
    [~,ix2] = min(abs(frametimes - gocue - tmax));

    frametimes = frametimes(ix1:ix2);

    %

    % traj
    feats = dat.obj.traj{view}(trial).featNames;
    feats2use = 4;
    jaw = normalize(dat.obj.traj{view}(trial).ts(ix1:ix2,2,feats2use),'range',[itrix itrix + 1]);

    feats2use = [1 2 3];
    tongue = squeeze(dat.obj.traj{view}(trial).ts(ix1:ix2,2,feats2use));
    % find all ix when tongue is visible
    tongvis = logical(sum(~isnan(tongue),2));

    jawtongvis = jaw;
    jawtongvis(~tongvis) = nan;
    
    plot(frametimes - gocue, jaw + dy,'Color',cols(2,:),'LineWidth',1.5);
    plot(frametimes - gocue, jawtongvis + dy,'k','LineWidth',1.5)
    
end
xlim([-1 1])
ylim([1 13])
xline(0,'k:')



