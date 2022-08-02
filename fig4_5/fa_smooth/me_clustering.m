clear,clc,close all
addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMALs/SESSIONs TO LOAD
meta(1).anm = 'JEB7';
meta(1).date = '2021-04-29';

meta(end+1).anm = 'JEB7';
meta(end).date = '2021-04-30';

meta(end+1).anm = 'JEB6';
meta(end).date = '2021-04-18';

meta(end+1).anm = 'JGR2';
meta(end).date = '2021-11-16';

meta(end+1).anm = 'JGR2';
meta(end).date = '2021-11-17';

meta(end+1).anm = 'JGR3';
meta(end).date = '2021-11-18';

meta(end+1).anm = 'EKH3';
meta(end).date = '2021-08-11';

% meta(end+1).anm = 'EKH3';
% meta(end).date = '2021-08-07';

meta(end+1).anm = 'EKH1';
meta(end).date = '2021-08-07';

meta = assignDataPath(meta);

% ------------------------------------------------------------------ %
% get some parameters that are common to all sessions/analyses
dfparams = getDefaultParams();

%% NEURAL ACTIVITY

% data objects contain
% 1) processed obj
% 2) processing params
% 3) factor analysis results
% 4) params used for factor analysis
% 5) meta data that we don't really need

[obj,fa,params,me] = loadProcessedData(meta);


%% cluster motion energy

six = 1; % session

temp = me(six).data';
temp(isnan(temp)) = 0;

times = [-1.7 -0.1];
[~,ix1] = min(abs(obj(six).time - times(1)));
[~,ix2] = min(abs(obj(six).time - times(2)));

clutemp = temp(:,ix1:ix2);


k = 3; % number of clusters;

% [idx,V_temp,D_temp] = spectralcluster(temp,k);

idx = kmeans(clutemp,k);

%%

figure;
for i = 1:k
    nexttile;
    trix = find(idx==i);
%     imagesc(obj(six).time(ix1:ix2),1:numel(trix),clutemp(trix,:));
    imagesc(obj(six).time,1:numel(trix),temp(trix,:));

    perf(i) = sum(obj(six).bp.hit(trix)) / numel(trix);

end

jaw = zeros(numel(obj(six).time),obj(six).bp.Ntrials);
for i = 1:obj(six).bp.Ntrials
    ts = obj(six).traj{2}(i).ts(:,2,2);
    tsinterp = interp1(obj(six).traj{1}(i).frameTimes-0.5-2.5, ts, obj(six).time); 
    jaw(:,i) = tsinterp;
end

figure;
for i = 1:k
    nexttile;
    trix = find(idx==i);
%     imagesc(obj(six).time(ix1:ix2),1:numel(trix),clutemp(trix,:));
    imagesc(obj(six).time,1:numel(trix),jaw(:,trix)');

end















