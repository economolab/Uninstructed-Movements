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

trix = [10,16,17,19,95,134];

vidfns = {'JEB15_2022-07-27_cam_0_date_2022_07_27_time_14_54_31_v001.avi',... % trial 10
          'JEB15_2022-07-27_cam_0_date_2022_07_27_time_14_55_14_v001.avi',... % trial 16
          'JEB15_2022-07-27_cam_0_date_2022_07_27_time_14_55_22_v001.avi',... % trial 17
          'JEB15_2022-07-27_cam_0_date_2022_07_27_time_14_55_38_v001.avi',... % trial 19
          'JEB15_2022-07-27_cam_0_date_2022_07_27_time_15_05_43_v001.avi',... % trial 95
          'JEB15_2022-07-27_cam_0_date_2022_07_27_time_15_11_02_v001.avi'};   % trial 134;

vidpth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\vidplot\vids\';

for itrix = 1:numel(trix)

    clear v hFig hAx hIm
    v = VideoReader(fullfile(vidpth,vidfns{itrix}));

%     sampleim = imread('tone.png');
%     delayim = imread('timer.png');
%     gocueim = imread('greendot.png');


    % get obj trajectories
    view = 1; % side cam
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

    % motion energy
    plotme = interp1(dat.obj.time,dat.me.data(:,trial),frametimes-gocue);


    % null/potent

    plotnull = squeeze(dat.rez.N_null(:,trial,:));
    plotnull = interp1(dat.obj.time,plotnull,frametimes-gocue);
    plotnull = mean(plotnull.^2,2);

    plotpotent = squeeze(dat.rez.N_potent(:,trial,:));
    plotpotent = interp1(dat.obj.time,plotpotent,frametimes-gocue);
    plotpotent = mean(plotpotent.^2,2);


    %

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
        'Position',[v.Width*0.1 v.Height*0.278 v.Width v.Height ],...
        'NextPlot','add',...
        'Visible','off',...
        'XTick',[],...
        'YTick',[]);
    hTraceAx = axes('Parent',hFig,...
        'Units','pixels',...
        'Position',[v.Width*0.08 0 292 v.Height*0.29 ],... % v.Width*1.2
        'NextPlot','add',...
        'Visible','on',...
        'XTick',[],...
        'YTick',[],...
        'Color','k');
    hIm = image(uint8(zeros(v.Height,v.Width,3)),...
        'Parent',hAx);

    % set up data
    msm = 4;
    npsm = 41;

    tempme = normalize(mySmooth(plotme,msm),'range',[5 35]);
    tempnull = normalize(mySmooth(plotnull,npsm),'range',[80 110]);
    temppotent = normalize(mySmooth(plotpotent,npsm),'range',[40 70]);

    lw = 1;
    hh(1) = plot(hTraceAx,frametimes,tempme,'y','LineWidth',1);
    hh(2) = plot(hTraceAx,frametimes,tempnull,'Color',[62, 168, 105]./255,'LineWidth',1);
    hh(3) = plot(hTraceAx,frametimes,temppotent,'Color',[255, 56, 140]./255,'LineWidth',1);
    hh(4) = xline(hTraceAx,gocue,'g-','LineWidth',2,'Alpha',0.7,'Color',[57, 205, 250]./255);

    % xlims = hTraceAx.XLim;
    text(hTraceAx,4,20,'ME','Color','y','FontSize',9,'FontWeight','bold');
    text(hTraceAx,4,55,'Potent','Color',[255, 56, 140]./255,'FontSize',9,'FontWeight','bold');
    text(hTraceAx,4,100,'Null','Color',[62, 168, 105]./255,'FontSize',9,'FontWeight','bold');


    ylims = hTraceAx.YLim;

    % loop through frames
    for iframe = 1:size(frames,4)
        im = frames(:,:,:,iframe);
        hIm.CData = flipud(im);

        % draw cursor
        traceline = xline(frametimes(iframe),'w','LineWidth',2);


        % epoch annotation
        t = '';
        if frametimes(iframe) >= sample && frametimes(iframe) < delay
            t = 'sample';
        elseif frametimes(iframe) >= delay && frametimes(iframe) < gocue
            delete(tt);
            t = 'delay';
        elseif frametimes(iframe) >= gocue && frametimes(iframe) <= (gocue + 0.25)
            t = 'go cue';
        end
        tt = text(hAx,15,220,t,'Color','w','FontSize',12);

        drawnow
        % Save the frame in structure for later saving to video file
        a = getframe(hAx);
        b = getframe(hTraceAx);
        s(iframe).cdata = cat(1,a.cdata,b.cdata);
        s(iframe).colormap = [];
        s(iframe).cdata(236:238,:,:) = 0;

        delete(tt);
        delete(traceline);
    end



    %
    % Remove any unused structure array elements

    % Write to a video file
    % This could be done within the original loop, but I wanted to show it
    % separately
    fnout = [anm '_' date '_trial' num2str(trial) '_alldim_ssm' '_100fps' '.mp4'];
    vOut = VideoWriter(fullfile('saved',fnout),'MPEG-4');
    vOut.FrameRate = 100;
    vOut.Quality = 100;
    open(vOut)
    for k = 1:numel(s)
        writeVideo(vOut,s(k))
    end
    close(vOut)


end




























