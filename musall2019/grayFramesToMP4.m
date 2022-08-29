
if numel(size(grayframes)) == 3
    temp = reshape(grayframes,size(grayframes,1),size(grayframes,2),nFrames,numel(params.trials2use));
else 
    temp = grayframes;
end

trix = randsample(size(temp,4),1);

figure;
for i = 1:size(temp,3)
    clf
    imagesc(temp(:,:,i,trix));
    drawnow
end


video = VideoWriter('exampleTrial_Cam1.avi'); %create the video object
video.FrameRate = 100;
open(video); %open the file for writing
for ii=1:size(temp,3) %where N is the number of images
  I = temp(:,:,ii,trix);
  writeVideo(video,I); %write the image to file
end
close(video); %close the file