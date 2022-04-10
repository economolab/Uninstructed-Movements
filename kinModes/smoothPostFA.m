function dat = smoothPostFA(dat,obj,params)

fs = 1./mean(diff(obj.time));                   % Find the sampling frequency for the time vector
% Define and implement a Butterworth filter for the single-trial data
[b,a] = butter(2,params.fcut_post_fa./fs./2);           % Butterworth filter: 2nd order, cutoff frequency determined using params.fcut and the sampling rate
% b = filter coefficient vector of numerator
% a = filter coefficient vector of denominators

dat.rates = filtfilt(b, a, dat.rates);                    % Digital filtering of the single-trial PSTHs with the filter described by vectors A and B

dat.factors = reshape(dat.factors,size(dat.rates,1),size(dat.rates,3),size(dat.factors,2));
dat.factors = permute(dat.factors,[1 3 2]);
dat.factors = filtfilt(b,a,dat.factors);

end