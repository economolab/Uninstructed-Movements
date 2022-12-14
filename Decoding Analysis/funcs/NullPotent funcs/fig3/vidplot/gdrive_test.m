

%%
gdrive = 'G:\My Drive\Economo-Lab\ExperimentsArchive\MAH14\Video\2022-08-15\Cam0';
vidfn = 'MAH14_2022-08-15_cam_0_date_2022_08_15_time_11_58_21_v001.avi'; 

clear v hFig hAx hIm
v = VideoReader(fullfile(gdrive,vidfn));    




nFrames = size(v.NumFrames,4);
s(nFrames) = struct('cdata',[],'colormap',[]);


% Set up figure, axes to hold image, image, and plot
% This allows us to keep everything positioned in the proper order, and
% just change the image data and the plot location every loop iteration
hFig = figure('MenuBar','none',...
    'Units','pixels','Color',[0 0 0],...
    'Position',[658 404 v.Width-10 v.Height-10]);

hAx = axes('Parent',hFig,...
    'Units','pixels',...
    'Position',[0 0 v.Width v.Height ],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);

hIm = image(uint8(zeros(v.Height,v.Width,3)),...
    'Parent',hAx);



for iframe = 1:v.NumFrames
    im = read(v,iframe);
%     im = frames(:,:,:,iframe);
    hIm.CData = flipud(im);

   
    drawnow
    % Save the frame in structure for later saving to video file
    s(iframe) = getframe(hAx);

end



%%
% Remove any unused structure array elements


% % Open a new figure and play the movie from the structure
% hFig2 = figure;
% movie(hFig2,s,1,v.FrameRate);

% Write to a video file
% This could be done within the original loop, but I wanted to show it
% separately
% fnout = [anm '_' date '_trial' num2str(trial) '_dim' num2str(dim) 'ssm_100fps' '.mp4'];
vOut = VideoWriter(fnout,'MPEG-4');
vOut.FrameRate = 100;
vOut.Quality = 100;
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k))
end
close(vOut)






























