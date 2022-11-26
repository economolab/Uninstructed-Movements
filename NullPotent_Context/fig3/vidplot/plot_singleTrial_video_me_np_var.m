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

% vidfn = 'JEB15_2022-07-27_cam_0_date_2022_07_27_time_15_05_43_v001.avi'; % trial 95
vidfn = 'JEB15_2022-07-27_cam_0_date_2022_07_27_time_14_55_14_v001.avi'; % trial 16

clear v hFig hAx hIm
v = VideoReader(vidfn);    

%%



% get obj trajectories
view = 1; % side cam
trial = 16;

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

%%

% motion energy
plotme = interp1(dat.obj.time,dat.me.data(:,trial),frametimes-gocue);


% null/potent
ndims = dat.rez.dPrep;
for i = 1:ndims
    temp = squeeze(dat.rez.N_null(:,trial,i));
    plotnull(:,i) = interp1(dat.obj.time,temp,frametimes-gocue);
end

ndims = dat.rez.dMove;
for i = 1:ndims
    temp = squeeze(dat.rez.N_potent(:,trial,i));
    plotpotent(:,i) = interp1(dat.obj.time,temp,frametimes-gocue);
end


% % jaw
% jawx = dat.obj.traj{view}(trial).ts(:,1,4);
% jawy = dat.obj.traj{view}(trial).ts(:,2,4);
% jawy = v.Height - jawy;
% plotme = jawy(ix1:ix2);

%%%%%
% figure; 
% subplot(2,1,1); hold on;
% plot(dat.obj.time,normalize(dat.me.data(:,trial),'range',[-2 6]))
% plot(dat.obj.time,dat.rez.N_null(:,trial,nullix))
% plot(dat.obj.time,dat.rez.N_potent(:,trial,potentix))
% subplot(2,1,2); hold on;
% plot(frametimes,normalize(plotme,'range',[-2 6]))
% plot(frametimes,plotnull)
% plot(frametimes,plotpotent)

%%

close all

nFrames = size(frames,4);
s(nFrames) = struct('cdata',[],'colormap',[]);


% Set up figure, axes to hold image, image, and plot
% This allows us to keep everything positioned in the proper order, and
% just change the image data and the plot location every loop iteration
hFig = figure('MenuBar','none',...
    'Units','pixels',...
    'Position',[100 100 v.Width v.Height]);
hAx = axes('Parent',hFig,...
    'Units','pixels',...
    'Position',[0 0 v.Width v.Height],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);
hIm = image(uint8(zeros(v.Height,v.Width,3)),...
    'Parent',hAx);
% hFeat1(1) = scatter(hAx,feat1x(1),feat1y(1),20,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
% hFeat2(1) = scatter(hAx,feat2x(1),feat2y(1),20,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
% hFeat3(1) = scatter(hAx,feat3x(1),feat3y(1),20,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
% numPrev = 50;
% for i = 1:numPrev
%     hFeat1(1+i) = scatter(hAx,feat1x(1),feat1y(1),20/(i+numPrev),'b','MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
%     hFeat2(1+i) = scatter(hAx,feat2x(1),feat2y(1),20/(i+numPrev),'b','MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
%     hFeat3(1+i) = scatter(hAx,feat3x(1),feat3y(1),20/(i+numPrev),'b','MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
% end


%

% create overlay

k = 1;
xmin = 0;
xmax = 1;
ms = 40;

nprev = 85;

msm = 7;
npsm = 11;

tempme = normalize(mySmooth(plotme,msm),'range',[5 35]);

meannull = var(plotnull,[],2);
tempnull = normalize(mySmooth(meannull,npsm),'range',[80 110]);

meanpotent = var(plotpotent,[],2);
temppotent = normalize(mySmooth(meanpotent,npsm),'range',[40 70]);

text(hAx,100,15,'Motion Energy','Color','y','FontSize',9);
text(hAx,100,50,'Potent','Color',[255, 56, 140]./255,'FontSize',9);
text(hAx,100,95,'Null','Color',[62, 168, 105]./255,'FontSize',9);

for iframe = 1:size(frames,4)
    im = frames(:,:,:,iframe);
    hIm.CData = flipud(im);

    % draw np
    if iframe > nprev
        hh(1) = plot(hAx,1:nprev,tempme(iframe-nprev+1:iframe),'y','LineWidth',2);
        hh(2) = plot(hAx,1:nprev,tempnull(iframe-nprev+1:iframe),'Color',[62, 168, 105]./255,'LineWidth',2);
        hh(3) = plot(hAx,1:nprev,temppotent(iframe-nprev+1:iframe),'Color',[255, 56, 140]./255,'LineWidth',2);
        
%         b = shadedErrorBar(1:nprev,meannull(iframe-nprev+1:iframe),errnull(iframe-nprev+1:iframe),{'Color',[62, 168, 105]./255,'LineWidth',1},0.5,hAx);
%         c = shadedErrorBar(1:nprev,meanpotent(iframe-nprev+1:iframe),errpotent(iframe-nprev+1:iframe),{'Color',[255, 56, 140]./255,'LineWidth',1},0.5,hAx);

    else
        hh(1) = plot(hAx,1:iframe+1,tempme(1:iframe+1),'y','LineWidth',2);
        hh(2) = plot(hAx,1:iframe,tempnull(1:iframe),'Color',[62, 168, 105]./255,'LineWidth',2);
        hh(3) = plot(hAx,1:iframe,temppotent(1:iframe),'Color',[255, 56, 140]./255,'LineWidth',2);

%         b = shadedErrorBar(1:iframe+1,meannull(1:iframe+1),errnull(1:iframe+1),{'Color',[62, 168, 105]./255,'LineWidth',1},0.5,hAx);
%         c = shadedErrorBar(1:iframe+1,meanpotent(1:iframe+1),errpotent(1:iframe+1),{'Color',[255, 56, 140]./255,'LineWidth',1},0.5,hAx);
    end


    % epoch annotation
    t = '';
    if frametimes(iframe) >= sample && frametimes(iframe) < delay
        t = 'sample';
    elseif frametimes(iframe) >= delay && frametimes(iframe) < gocue
        t = 'delay';
    elseif frametimes(iframe) >= gocue && frametimes(iframe) <= (gocue + 0.25)
        t = 'go cue';
    end
    tt = text(hAx,15,220,t,'Color','w','FontSize',12);

    
    drawnow
    % Save the frame in structure for later saving to video file
    s(iframe) = getframe(hAx);

    delete(tt);
    delete(hh);
%     delete(a); 
%     delete(b.mainLine); delete(b.patch); delete(b.edge);
%     delete(c.mainLine); delete(c.patch); delete(c.edge);
end



%%
% Remove any unused structure array elements


% % Open a new figure and play the movie from the structure
% hFig2 = figure;
% movie(hFig2,s,1,v.FrameRate);

% Write to a video file
% This could be done within the original loop, but I wanted to show it
% separately
fnout = [anm '_' date '_trial' num2str(trial) '_dim' num2str(dim) '_var_100fps' '.mp4'];
vOut = VideoWriter(fnout,'MPEG-4');
vOut.FrameRate = 100;
vOut.Quality = 100;
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k))
end
close(vOut)































