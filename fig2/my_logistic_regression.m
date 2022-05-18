clear,clc,close all
% IMPORT DATA
load fisheriris

% DATA SPECIFICATION
Y = onehotencode(categorical(species),2); % as many columns as there are unique categories, 1s corresponding to what category each row of species is.
Ynew = nan(size(Y,1),1);
for i = 1:size(Y,2)
    ix = Y(:,i) > 0;
    Ynew(ix) = Y(ix,i) * (i);
end
Y = Ynew;

X = meas;

% VARIABLE NAMES
Vnames = {'x1','x2','x3','x4'}; % covariate/regressor names, for fisheriris these are like pedal length,wdith, and sepal length/width

% STARTING VALUES
b = zeros(length(Vnames),1); % weights

% OPTIMISATION
options = optimoptions(@fminunc,'Display','iter','MaxIterations',1e3,'MaxFunctionEvaluations',1e5);
tic;
[paramhat,fval,~,~,grad,hessian] = fminunc(@(b)loglike_logisticreg(b, X, Y), b, options);