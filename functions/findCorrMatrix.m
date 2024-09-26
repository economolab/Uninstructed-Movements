function corrmatrix = findCorrMatrix(modes)
nWindows = size(modes,2);
corrmatrix = NaN(nWindows,nWindows);
for tt = 1:nWindows                 % For every time point where you calculated a mode
    currmode1 = modes(:,tt);        % Get the mode found for the current time point
    for ww = 1:nWindows
        currmode2 = modes(:,ww);    % Get the modes found at all other time points
        R = corrcoef(currmode1,currmode2);      % Find the correlation between the first mode and all of the other modes found
        corrmatrix(tt,ww) = R(1,2);
    end
end
end