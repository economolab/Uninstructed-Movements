%% INPUTS: 
% selectNorm = matrix of selectivity values for all cells across trials
% sort_by = time period that you want to sort the selectivity heatmap by

% order = 'first' indicates that this is the first heatmap you will be
% sorting and plotting.  So you want to save the order in which the cells
% are sorted...'second' indicates that this is not the first heatmap you
% are sorting and that you want to sort by a previously determined order

% sIx = the indices by which you want to sort; omit this if order = 'first'

function [selectNorm,sIx] = sortSelectivity(selectNorm,sort_by,order,sIx)

if strcmp(order,'first')
    [selectSorted,sIx] = sort(selectNorm(sort_by,:),2,'descend'); selectNorm = selectNorm(:,sIx);  % Sort by the selectivity at specified time point
else
    selectNorm = selectNorm(:,sIx);

end