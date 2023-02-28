% Quantifying behavioral performance and cortical dependence in the
% Alternating Context Task
% -------------------------------------------------------------------------------------
clear,clc,close all

% add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 1';
addpath(genpath(fullfile(figpth,'Utils')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% set conditions to calculate behavioral performance for 
params.condition(1) = {'L&~stim.enable&~autowater'};              % right, no stim, 2AFC
params.condition(end+1) = {'R&~stim.enable&~autowater'};              % left, no stim, 2AFC
params.condition(end+1) = {'L&stim.enable&~autowater'};               % right, stim, afc
params.condition(end+1) = {'R&stim.enable&~autowater'};               % left, stim, afc

params.condition(end+1) = {'L&~stim.enable&autowater'};              % right, no stim, aw
params.condition(end+1) = {'R&~stim.enable&autowater'};          % left, no stim, aw
params.condition(end+1) = {'L&stim.enable&autowater'};               % right, stim, aw
params.condition(end+1) = {'R&stim.enable&autowater'};               % left, stim, aw
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];
meta = loadMAH13_MCStim(meta,datapth);
meta = loadMAH14_MCStim(meta,datapth);

%% LOAD DATA

% ----------------------------------------------
% -- Behavioral and Video Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadBehavSessionData(meta,params);
%%
close all

sessix = 2; 

cols = getColors_Updated();
clrs.rhit = cols.rhit;
clrs.lhit = cols.lhit;

conds = [1 2 3 4];

figure(1)
plotLickRaster(sessix,clrs,obj,params,'2AFC',conds);

%%
clrs.rhit = cols.rhit_aw;
clrs.lhit = cols.lhit_aw;

conds = [5 6 7 8];

figure(2)
plotLickRaster(sessix,clrs,obj,params,'AW',conds);