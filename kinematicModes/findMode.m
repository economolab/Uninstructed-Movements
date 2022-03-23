function [mode, dat] = findMode(obj, feat, params)

fr = cat(3, obj.trialpsth_cond{params.cond});

fs = 1./mean(diff(obj.time));
[b,a] = butter(2,params.fcut./fs./2);
filtfr = filtfilt(b, a, fr);

filtfr = permute(filtfr(params.tix,:, :), [ 1 3 2]);
filtfr = reshape(filtfr, size(filtfr, 1)*size(filtfr, 2), size(filtfr, 3));

if params.fa
[~,~,~,~,dat]  = factoran(filtfr,10);
else
    dat = filtfr;
end

% code for finding mode by xcorr

feat = feat(params.tix, :);
feat = feat(:);
feat(isnan(feat)) = 0;

switch params.method
    case 'xcorr'
        mode = doXCorr(dat, feat);
       
    case 'regress'
        mode = regress(feat,dat);
       
        mode = mode./sum(abs(mode));
    otherwise
        mode = zeros(1, 1);
        return;
        
end

dat = reshape(dat, numel(params.tix), size(fr, 3), size(dat, 2));
dat = permute(dat, [1 3 2]);

