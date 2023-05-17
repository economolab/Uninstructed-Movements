
clear,clc,close all
%% data object and video file

datapth = 'C:\Users\Jackie Birnbaum\Documents\Data\DataObjects\JEB13';         % Directory where data obj is stored
objfn = 'data_structure_JEB13_2022-09-13';                                     % File name of data object
vidfn = 'JEB13_2022-09-13_cam_0_date_2022_09_13_time_09_34_56_v001.avi';       % Filename of video
load(fullfile(datapth, objfn));        % Load the data object

v = VideoReader(fullfile(datapth, vidfn));             

%%
% get obj trajectories
view = 1; % side cam
feat = 3; % nose
trial = 16;

% nose
feat1x = obj.traj{view}(trial).ts(:,1,3); % trial 24 corresp0onds to video file
feat1y = obj.traj{view}(trial).ts(:,2,3); % trial 24 corresp0onds to video file
feat1y = v.Height - feat1y; % y pos starts from top of image

% tongue
feat2x = obj.traj{view}(trial).ts(:,1,1);
feat2y = obj.traj{view}(trial).ts(:,2,1);
feat2y = v.Height - feat2y;
feat2x(1:1000) = nan;
feat2y(1:1000) = nan;

% jaw
feat3x = obj.traj{view}(trial).ts(:,1,2);
feat3y = obj.traj{view}(trial).ts(:,2,2);
feat3y = v.Height - feat3y;

%%
% Preallocate structure to store video frames
% Note that this is just an initial guess for preallocation based on
% duration and framerate, the video may have fewer or more frames
nFrames = ceil(v.FrameRate*v.Duration);
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
hFeat1(1) = scatter(hAx,feat1x(1),feat1y(1),20,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
hFeat2(1) = scatter(hAx,feat2x(1),feat2y(1),20,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
hFeat3(1) = scatter(hAx,feat3x(1),feat3y(1),20,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
numPrev = 10;
for i = 1:numPrev
    hFeat1(1+i) = scatter(hAx,feat1x(1),feat1y(1),20/(i+numPrev),'b','MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
    hFeat2(1+i) = scatter(hAx,feat2x(1),feat2y(1),20/(i+numPrev),'b','MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
    hFeat3(1+i) = scatter(hAx,feat3x(1),feat3y(1),20/(i+numPrev),'b','MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
end
% drawing the trace works, but need to delete the last hDot at each
% iteration

%%
% Loop through video, grabbing frames and updating plots

% set up colors
stretch = 0.1;
cols1 = (cool(ceil(nFrames*stretch)));
temp = cols1;
for i = 1:(1/stretch)
    if mod(i,2)~=0
        temp = [temp ; flipud(cols1)];
    else
        temp = [temp ; cols1];
    end
end
cols1 = temp;

cols2 = (spring(ceil(nFrames*stretch)));
temp = cols2;
for i = 1:(1/stretch)
    if mod(i,2)~=0
        temp = [temp ; flipud(cols2)];
    else
        temp = [temp ; cols2];
    end
end
cols2 = temp;

cols3 = (summer(ceil(nFrames*stretch)));
temp = cols3;
for i = 1:(1/stretch)
    if mod(i,2)~=0
        temp = [temp ; flipud(cols3)];
    else
        temp = [temp ; cols3];
    end
end
cols3 = temp;

%% create overlay

k = 1;
init = v.CurrentTime;
xmin = 0;
xmax = 1;
ms = 40;
while hasFrame(v)
    im = readFrame(v);
    hIm.CData = flipud(im);
    % draw annotation
    hFeat1(1) = scatter(hAx,feat1x(k),feat1y(k),ms,'MarkerFaceColor',cols1(k,:),'MarkerEdgeColor',cols1(k,:));
    hFeat2(1) = scatter(hAx,feat2x(k),feat2y(k),ms,'MarkerFaceColor',cols2(k,:),'MarkerEdgeColor',cols2(k,:));
    hFeat3(1) = scatter(hAx,feat3x(k),feat3y(k),ms,'MarkerFaceColor',cols3(k,:),'MarkerEdgeColor',cols3(k,:));
    % draw previous annotations
    if k>numPrev
        for i = 1:numPrev
            hFeat1(i+1) = scatter(hAx,feat1x(k-i),feat1y(k-i),ceil(ms/(i-0.5)),'MarkerFaceColor',cols1(k-i,:), ...
                'MarkerEdgeColor',cols1(k-i,:),'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
            hFeat2(i+1) = scatter(hAx,feat2x(k-i),feat2y(k-i),ceil(ms/(i-0.5)),'MarkerFaceColor',cols2(k-i,:), ...
                'MarkerEdgeColor',cols2(k-i,:),'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
            hFeat3(i+1) = scatter(hAx,feat3x(k-i),feat3y(k-i),ceil(ms/(i-0.5)),'MarkerFaceColor',cols3(k-i,:), ...
                'MarkerEdgeColor',cols3(k-i,:),'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
        end
    end
    
    drawnow
    % Save the frame in structure for later saving to video file
    s(k) = getframe(hAx);
    k = k+1;
    delete(hFeat1);
    delete(hFeat2);
    delete(hFeat3);
end
v.CurrentTime = init; % reset time to play back from first frame

%%
% Remove any unused structure array elements
s(k:end) = [];

% % Open a new figure and play the movie from the structure
% hFig2 = figure;
% movie(hFig2,s,1,v.FrameRate);

% Write to a video file
% This could be done within the original loop, but I wanted to show it
% separately
vOut = VideoWriter('JGR2_trial24_tongue_jaw_nose_100.mp4','MPEG-4');
vOut.FrameRate = 100;
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k))
end
close(vOut)









