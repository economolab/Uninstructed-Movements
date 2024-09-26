% Script for adjusting spike times in data object according to the buffer
% time from SpikeGLX
clear; clc;

load('C:\Users\Jackie Birnbaum\Documents\Data\DataObjects\JEB19\data_structure_JEB19_2023-04-21.mat')

nClu = size(obj.clu{1},2);
nTrials = obj.bp.Ntrials;
probenum = 1;
for clu = 1:nClu
    spktms = obj.clu{probenum}(clu).trialtm;
    spktms = spktms-0.5;
    obj.clu{probenum}(clu).trialtm = spktms;
end