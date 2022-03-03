% INPUTS: 
% taxis = the time axis that you want to be on the xaxis
% sort_by = the time point in the trial that you want to sort selectivity
% at (i.e. lateDelay or alignEvent)
% numCells_tot = the number of cells that you are plotting selectivity for
% selectNorm = selectivity for each cell (can be direction or context
% selectivity), normalized to the maximum selectivity value found in any of
% the cells)

% OUTPUTS:
% A heatmap of selectivity (ranges from -1 to 1), aligned to the align
% event

function plotSelectivityHeatmap(taxis,sort_by,numCells_tot,selectNorm)
[selectSorted,sIx] = sort(selectNorm(sort_by,:),2,'descend'); selectNorm = selectNorm(:,sIx);  % Sort by the selectivity at specified time point

alignEv = find(taxis == 0);       
earlyDelay = find(taxis>-0.855 & taxis<-0.845);     % t = 0.85 s before first lick
lateDelay = find(taxis>-0.253 & taxis<-0.247);      % t = 0.3 s before first lick
preSample = find(taxis<-2.495);                     % t = 2.5 s before first lick

if sort_by==alignEv                                 % Specify where to draw dotted line (at time that you are sorting by)
    l=0;
elseif sort_by==earlyDelay
    l=-0.85;
elseif sort_by==lateDelay
    l=-0.25;
elseif sort_by==preSample
    l=-2.5;
end
imagesc(taxis,1:numCells_tot,selectNorm');
hold on;
line([l,l],[1,numCells_tot],'Color','black','LineStyle','--')
ax = gca;
c = colorbar(ax);
c.FontSize = 14;
lim=[-1 1];
caxis(lim);
colormap(ax,flipud(jet))
ylabel(c,'L <---> R')
end