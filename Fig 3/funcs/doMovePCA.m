function [coeff, score] = doMovePCA(trialmove, Ndims)
% concatenate data into a 2D array of size (trials*Time,F) where F is number of
% features, T is time points
A = [];                                         % Pre-allocate activity matrix of move features
for i = 1:size(trialmove,2)                     % For all trials...
    A = [A; trialmove(:,i,:)];                  % Concatenate the activity matrices for all features (aka [tonguevel; jawvel; pawvel; ME])
end
A = squeeze(A);
A = fillmissing(A,"constant",0);
meancenterA = A-mean(A,'omitnan');
meancenterA = fillmissing(meancenterA,'constant',0);
meancenterA(isinf(meancenterA)) = 0;

% coeff is the principal components. Should be of size (Nfeatures,Ndims). Each column
% is a principal component, and each element of each column is the weight
% associated with a feature.
% score is the projection of the full dimensional move data onto the principal
% components (score = psth * coeff)
[coeff, score] = pca(meancenterA, 'NumComponents', Ndims);          % Perform PCA on centered data
end