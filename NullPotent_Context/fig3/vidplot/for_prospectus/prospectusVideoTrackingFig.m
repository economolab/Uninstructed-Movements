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

trix = 16;

% vidfns = {'JEB15_2022-07-27_cam_0_date_2022_07_27_time_14_55_14_v001.avi'}; % trial 16 side cam
vidfns = {'JEB15_2022-07-27_cam_1_date_2022_07_27_time_14_55_14_v001.avi'}; % trial 16 side cam
         
vidpth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\vidplot\vids\';

for itrix = 1:numel(trix)

    clear v hFig hAx hIm
    v = VideoReader(fullfile(vidpth,vidfns{itrix}));

%     sampleim = imread('tone.png');
%     delayim = imread('timer.png');
%     gocueim = imread('greendot.png');


    % get obj trajectories
    view = 2; % side cam
    trial = trix(itrix);

    frametimes = dat.obj.traj{view}(trial).frameTimes - 0.5;
    sample = dat.obj.bp.ev.sample(trial);
    delay = dat.obj.bp.ev.delay(trial);
    gocue = dat.obj.bp.ev.goCue(trial);


    % resmaple obj.time to frametimes, but first need to get frametimes from
    % -2.5 from go cue to 2.5 to go cue
    tmin = dat.params.tmin;
    tmax = dat.params.tmax;

    [~,ix1] = min(abs(frametimes - gocue - tmin));
    [~,ix2] = min(abs(frametimes - gocue - tmax));

    frames = read(v,[ix1 ix2]);
    frametimes = frametimes(ix1:ix2);

    %

    % traj
    feats = dat.obj.traj{view}(trial).featNames;
    if view == 2
        feats2use = [1 5 6 8 9 10];
    else
        feats2use = [2 4 6];
    end
    traj = dat.obj.traj{view}(trial).ts(ix1:ix2,1:2,feats2use);

    close all

    nFrames = size(frames,4);

    sm = 21;
    lw = 1;

    % loop through frames
    f = figure;
    frameix = 1101;
    cols = linspecer(numel(feats2use)+5);
    for iframe = frameix
        im = frames(:,:,:,iframe);
        imshow(im);
        title([anm '_' date '  Trial: ' num2str(trial) ' Frame: ' num2str(iframe)],'Interpreter','none')
        hold(gca,'on');
        for ifeat = 1:numel(feats2use)
            if view == 2
                colix = ifeat + 3;
            else
                colix = ifeat;
            end
            if strcmpi(feats{feats2use(ifeat)},'jaw')
                plot(traj(iframe,1,ifeat), traj(iframe,2,ifeat),'.','MarkerSize',35,'Color',cols(colix,:))
            else
                plot(traj(iframe,1,ifeat), traj(iframe,2,ifeat),'.','MarkerSize',35,'Color',cols(colix,:))
            end
        end

        hold(gca,'off')
        drawnow
%         pause
    end


    
    cols = linspecer(numel(feats2use));
    dy = 0;
    for ifeat = 1:numel(feats2use)
        f = figure; hold on;
        if strcmpi(feats{feats2use(ifeat)},'tongue')
            ix = ifeat + 1;
        else
            ix = ifeat;
        end
        plot(frametimes - gocue, traj(:,2,ix) + dy,'LineWidth',2,'Color',cols(ifeat,:))
        xline(0,'k:','LineWidth',1)
        xline(sample-gocue,'k:','LineWidth',1)
        xline(delay-gocue,'k:','LineWidth',1)
        ylim([60 215])
        xlim([-0.5 1.5])
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
    end



%     %
%     % Remove any unused structure array elements
% 
%     % Write to a video file
%     % This could be done within the original loop, but I wanted to show it
%     % separately
%     fnout = [anm '_' date '_trial' num2str(trial) '_alldim_ssm' '_100fps' '.mp4'];
%     vOut = VideoWriter(fullfile('saved',fnout),'MPEG-4');
%     vOut.FrameRate = 100;
%     vOut.Quality = 100;
%     open(vOut)
%     for k = 1:numel(s)
%         writeVideo(vOut,s(k))
%     end
%     close(vOut)


end




























