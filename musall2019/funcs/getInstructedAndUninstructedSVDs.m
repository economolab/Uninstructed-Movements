function [instructed, uninstructed] = getInstructedAndUninstructedSVDs(dat,ix,toPlot)

% % video SVD


temp = squeeze(mean(dat,2));
instructed = temp;

% instructed_vidR(1:goCueIx-1,:) = 0; % trial-averaged video svd, everything before go cue zeroed out
instructed = repmat(instructed,1,1,size(dat,2)); % make nTrials copies in 3rd dimension, adding trials back
instructed = permute(instructed,[1 3 2]);

if toPlot
    figure;
    subplot(1,2,1)
    imagesc(obj.time,1:nrDims,temp')
    hold on
    xline(obj.time(ix),'k--','LineWidth',2)
    subplot(1,2,2)
    imagesc(obj.time,1:nrDims,instructed')
    hold on
    xline(obj.time(ix),'k--','LineWidth',2)
end


% subtract instructed vid from vidR to get uninstructed component
% uninstructed_vidR = vidR - instructed_vidR;
uninstructed = dat;

end