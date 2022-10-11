function [xvel, yvel] = findVelocity(xpos, ypos, feat)

for i = 1:size(xpos,2) % for each trial
    tsinterp = [xpos(:,i) ypos(:,i)];
    basederiv = median(diff(tsinterp),'omitnan');

    % find the difference between the feat velocity and the
    % baseline feature velocity (NOT FOR TONGUE)
    xvel(:,i) = gradient(tsinterp(:,1));
    yvel(:,i) = gradient(tsinterp(:,2));
    if ~contains(feat,'tongue') %|| strcmp(traj(trix).featNames{feat},'nose')
        xvel(:,i) = xvel(:,i) - basederiv(1);
        yvel(:,i) = yvel(:,i) - basederiv(1);
    end

    % fill missing values for all features except tongue
    if ~contains(feat,'tongue')
        xvel(:,i) = fillmissing(xvel(:,i),'nearest');
        yvel(:,i) = fillmissing(yvel(:,i),'nearest');
    end
end


end  % findVelocity





