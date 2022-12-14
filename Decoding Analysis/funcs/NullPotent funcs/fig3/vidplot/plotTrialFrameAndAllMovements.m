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

trix = [10];

view = 1;

if view == 1
    vidfns = {'JEB15_2022-07-27_cam_0_date_2022_07_27_time_14_54_31_v001.avi'}; % trial 10
else
    vidfns = {'JEB15_2022-07-27_cam_1_date_2022_07_27_time_14_54_31_v001.avi'}; % trial 10
end
        

vidpth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\vidplot\vids\';

for itrix = 1:numel(trix)

    clear v hFig hAx hIm traj
    v = VideoReader(fullfile(vidpth,vidfns{itrix}));

%     sampleim = imread('tone.png');
%     delayim = imread('timer.png');
%     gocueim = imread('greendot.png');


    % get obj trajectories
    trial = trix(itrix);

    frametimes = dat.obj.traj{view}(trial).frameTimes - 0.5;
    sample = dat.obj.bp.ev.sample(trial);
    delay = dat.obj.bp.ev.delay(trial);
    gocue = dat.obj.bp.ev.goCue(trial);


    % resmaple obj.time to frametimes, but first need to get frametimes from
    % -2.5 from go cue to 2.5 to go cue
    tmin = dat.params.tmin;
    tmax = -0.02;

    [~,ix1] = min(abs(frametimes - gocue - tmin));
    [~,ix2] = min(abs(frametimes - gocue - tmax));

    frames = read(v,[ix1 ix2]);
    frametimes = frametimes(ix1:ix2);


    % trajectories across all trials, before go cue
    if view == 1
        feats = [4];
    else
        feats = [5 6 8 9 10];
    end
    for iii = 1:dat.obj.bp.Ntrials
        traj(:,:,:,iii) = dat.obj.traj{view}(iii).ts(ix1:ix2,[1 2],feats);
    end
    traj = permute(traj, [1 4 2 3]);
    traj = reshape(traj,size(traj,1)*size(traj,2),size(traj,3),size(traj,4));
    

    close all

    nFrames = size(frames,4);
    s(nFrames) = struct('cdata',[],'colormap',[]);


    % Set up figure, axes to hold image, image, and plot
    % This allows us to keep everything positioned in the proper order, and
    % just change the image data and the plot location every loop iteration
    hFig = figure('MenuBar','none',...
        'Units','pixels','Color',[0 0 0],...
        'Position',[658 404 v.Width*1.15 v.Height*1.25]);
    

    hAx = axes('Parent',hFig,...
        'Units','pixels',...
        'Position',[v.Width*0.05 v.Height*0.05 v.Width v.Height ],...
        'NextPlot','add',...
        'Visible','off',...
        'XTick',[],...
        'YTick',[]);
%     hAx = axes('Parent',hFig,...
%         'Units','pixels',...
%         'Position',[v.Width v.Height v.Width v.Height ],...
%         'NextPlot','add',...
%         'Visible','off',...
%         'XTick',[],...
%         'YTick',[]);
    hIm = image(uint8(zeros(v.Height,v.Width,3)),...
        'Parent',hAx);

    cols = linspecer(3,'qualitative');


    iframe = 192;
    im = frames(:,:,:,iframe);
    hIm.CData = flipud(im);

    for fix = 1:numel(feats)
        xs = round(squeeze(traj(:,1,fix))) + 15;
        ys = v.Height - round(squeeze(traj(:,2,fix)));

        ff = plot(hAx, xs, ys , 'Color', cols(fix+1,:),'LineWidth',1);
        ff.Color = [ff.Color 0.2];
    end

    
    fnout = [anm '_' date '_trial' num2str(trial) 'view_' num2str(view) '_prospectus'];
    saveas(hFig,fnout,'svg')




end




























