%% Inseparable case

% Latent dynamics (L1)
close all
xx = 0:1/1e4:0.5;
L1 = square(2*pi*10*xx,40);
L1(L1==-1) = 0;

onIX = (L1==1);
start = []; 
stop = 1;
for ii = 1:(length(onIX)-1)
    if onIX(ii) && ~onIX(ii+1)
        start = [start,ii];
    elseif ~onIX(ii) && onIX(ii+1)
        stop = [stop,ii];
    end
end


cog = L1;
motIn = L1;
for ss = 1:length(start)
    ix1 = stop(ss);
    ix2 = start(ss);

    % Cognitive process (C)
    dur = linspace(0,1,length(ix1:ix2));
    cog(ix1:ix2) = dur;

    % Motor process (M)
    Fr = 25;        % Change this param to alter frequency of sine wave in "on" periods of motor process
    len = linspace(0,Fr,length(dur));
    motIn(ix1:ix2) = sin(len);
end

subplot(3,1,1)
plot(xx,L1+9); hold on;
plot(xx,cog+7)
plot(xx,motIn+5)
ylim([-1 11])
legend('L1','C','M')
title('Inseparable')

clearvars -except xx cog mot L1
%% Separable and independent

% Generate L2
L2 = zeros(1,length(L1));
startPct = 0.15;            % How far into the trace you want L2 to "turn on"
stopPct = 0.75;             % How far into the trace you want L2 to "turn off"
startIX = floor(startPct*length(L2));
stopIX = floor(stopPct*length(L2));
L2(startIX:stopIX) = 1;

clearvars -except xx cog mot L1 L2

motInd = L2;
onIX = find(L2==1);
startIX = onIX(1);
stopIX = onIX(end);
dur = linspace(0,1,length(startIX:stopIX));

% Motor process (M)
Fr = 103.5;        % Change this param to alter frequency of sine wave in "on" periods of motor process
len = linspace(0,Fr,length(dur));
motInd(startIX:stopIX) = sin(len);


subplot(3,1,2)
plot(xx,L1+9); hold on;
plot(xx,cog+7)
plot(xx,L2+5)
plot(xx,motInd+3)
ylim([-1 11])
legend('L1','C','L2','M')
title('Separable, independent')

clearvars -except xx cog mot L1 L2 motInd

%% Separable and dependent

L2dep = L2;

onIX = find(L2==1);
startIX = onIX(1);
stopIX = onIX(end);
dur = linspace(0,1,length(startIX:stopIX));

motDep = motInd;
Fr = 103.5;        % Change this param to alter frequency of sine wave in "on" periods of motor process
len = linspace(0,Fr,length(dur));
motDep(startIX:stopIX) = 0.5*sin(len);
clear dur

onL1 = (L1==1);
start = []; 
stop = 1;
for ii = 1:(length(onL1)-1)
    if onL1(ii) && ~onL1(ii+1)
        start = [start,ii];
    elseif ~onL1(ii) && onL1(ii+1)
        stop = [stop,ii];
    end
end


for ss = 1:length(start)
    ix1 = stop(ss);
    ix2 = start(ss);

    % Generate L2
    L2dep(ix1:ix2) = 0.5+L2dep(ix1:ix2);

    % Generate Motor process (M)
    dur = ix1:ix2;
    Fr = 24;        % Change this param to alter frequency of sine wave in "on" periods of motor process
    len = linspace(0,Fr,length(dur));
    motDep(ix1:ix2) = sin(len);
    
end

L2dep(L2==0) = 0;
motDep(L2==0) = 0;


% motInd = L2;
% onIX = find(L2==1);
% startIX = onIX(1);
% stopIX = onIX(end);
% dur = linspace(0,1,length(startIX:stopIX));
% 
% % Motor process (M)
% Fr = 104.5;        % Change this param to alter frequency of sine wave in "on" periods of motor process
% len = linspace(0,Fr,length(dur));
% motInd(startIX:stopIX) = sin(len);


subplot(3,1,3)
plot(xx,L1+9); hold on;
plot(xx,cog+7)
plot(xx,L2dep+5)
plot(xx,motDep+3)
ylim([-1 11])
legend('L1','C','L2','M')
title('Separable, dependent')

