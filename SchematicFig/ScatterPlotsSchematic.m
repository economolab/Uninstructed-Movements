%%% Code for generating schematic plots
xmax = 300;
cmax = xmax/2;
buffer = 0;
cog1 = buffer:(cmax-buffer/2);
cog2 = (cmax+buffer/2):(xmax-buffer);
%% Inseparable
yval1 = 80;
noiselevel = 3;
mot1 = yval1*ones(1,length(cog1));
temp = linspace((-1*noiselevel),noiselevel,length(mot1));
noise = randsample(temp,length(mot1));
motnoise1 = mot1+noise;

yval2 = 220;
mot2 = yval2*ones(1,length(cog2));
noise = randsample(temp,length(mot1));
motnoise2 = mot2+noise;

axlim = [45 255];

figure()
scatter(cog1,motnoise1,20,'black','filled'); hold on;
scatter(cog2,motnoise2,20,'black','filled')
xlim(axlim)
ylim(axlim)
ax = gca;
ax.TickDir = 'out';
xlabel('Cognitive neural dimension')
ylabel('Motor neural dimension')
title('Inseparable')
%% Separable and independent
temp = linspace(50,250,length(mot1));
mot1 = randsample(temp,length(mot1));

mot2 = randsample(temp,length(mot2));

figure()
scatter(cog1,mot1,20,'black','filled'); hold on;
scatter(cog2,mot2,20,'black','filled')
xlim(axlim)
ylim(axlim)
ax = gca;
ax.TickDir = 'out';
xlabel('Cognitive neural dimension')
ylabel('Motor neural dimension')
title('Separable and independent')
%%
clearvars -except xmax cmax buffer cog1 cog2 axlim
%% Separable but dependent
yval1 = 80;
dependlevel = 0.6;
dependmots = floor(dependlevel*(length(cog1)));
randmots = length(cog1)-dependmots;
noiselevel = 9;

emb1 = yval1*ones(1,dependmots);
temp = linspace((-1*noiselevel),noiselevel, length(emb1));
dep = randsample(temp,dependmots);
dep1 = emb1+dep;
emb1 = yval1*ones(1,randmots);
temp = linspace(-80,175,length(emb1));
noise = randsample(temp,randmots);
noise1 = emb1+noise;
mot1 = [dep1, noise1]; mot1 = randsample(mot1,length(mot1));

yval2 = 220;
emb2 = yval2*ones(1,dependmots);
temp = linspace((-1*noiselevel),noiselevel,length(emb2));
dep = randsample(temp,dependmots);
dep2 = emb2+dep;
emb2 = yval1*ones(1,randmots);
temp = linspace(-80,175,length(emb2));
noise = randsample(temp,randmots);
noise2 = emb2+noise;
mot2 = [dep2, noise2]; mot2 = randsample(mot2,length(mot2));


axlim = [45 255];

figure()
scatter(cog1,mot1,20,'black','filled'); hold on;
scatter(cog2,mot2,20,'black','filled')
xlim(axlim)
ylim(axlim)
ax = gca;
ax.TickDir = 'out';
xlabel('Cognitive neural dimension')
ylabel('Motor neural dimension')
title('Separable but dependent')