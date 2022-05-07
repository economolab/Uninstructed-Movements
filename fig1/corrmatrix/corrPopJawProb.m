clear,clc,close all

% correlation at each time point of the jaw prob corr matrix and neural pop
% corr matrix (correlate slices of the corr matrices)

dat = load('popcorrmatrix.mat');
pop = dat.corr_matrix_selectivity;

dat = load('jawprobcorrmatrix.mat');
jaw = dat.corr_matrix_selectivity;

%%

pop(isnan(pop)) = 0;
jaw(isnan(jaw)) = 0;

corrs = zeros(size(pop,1),1);
for i = 1:numel(corrs)
    temp = corrcoef(pop(:,i),jaw(:,i));
    corrs(i) = temp(1,2);
end

%%

time = -2.5:0.005:2.5;
time = time(1:end-1);
figure;
plot(time,corrs,'k','LineWidth',3.5)
xlabel('Time (s) from go cue')
ylabel('Correlation b/w jawprob and neural pop corr matrices')
ax = gca;
ax.FontSize = 20;
















