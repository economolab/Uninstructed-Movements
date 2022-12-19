% Function for finding kinematic modes using cross-correlation
% INPUTS:
% X = dat (time*trials x cells) --> predictors.  Each cell is a predictor
% that you are trying to predict the kinematic feature with
% Y = feat (time*trials x 1) --> observations. The kinematic feature you are trying
% to predict

% OUTPUTS:
% mode = (cells x 1) --> loadings/weight of each cell in being able to
% predict the feature

% Y are observations (N x 1)
% X are predictors (N x M)
% mode is loadings (M x 1)

function mode = doXCorr(X, Y)
r = zeros(size(X, 2), 1);       % (cells x 1)

for i = 1:size(X, 2)            % For every cell...
    x = X(:, i);                    % Get the PSTH for all times and trials for that cell
    tmp = corrcoef(x, Y);           % Find the correlation coefficient between the cell's firing rate and the feature measure
    r(i) = tmp(1,2);                % Store the corr coefficient
end

mode = r./sum(abs(r));          % Normalize all correlation coefficients to the sum of the absolute value of all coefficients