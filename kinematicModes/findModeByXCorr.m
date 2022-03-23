function mode = findModeByXCorr(obj, cond, feat, ix, fcut, method)

fr = cat(3, obj.trialpsth_cond{cond});

fs = 1./mean(diff(obj.time));
[b,a] = butter(2,fcut./fs./2);
filtfr = filtfilt(b, a, fr);

% code for finding mode by xcorr
filtfr = permute(filtfr(ix,:, :), [ 1 3 2]);
filtfr = reshape(filtfr, size(filtfr, 1)*size(filtfr, 2), size(filtfr, 3));

feat = feat(ix, :);
feat = feat(:);
feat(isnan(feat)) = 0;

switch method
    case 'xcorr'
        mode = doXCorr(filtfr, feat);
    case 'regression'
        
    otherwise
        mode = zeros(1, 1);
        return;
        
end

