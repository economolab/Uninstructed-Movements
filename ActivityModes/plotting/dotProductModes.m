function dotProductModes(rez,modes,titleStr)

[modenames,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
for i = 1:numel(modenames)
    modenames{i} = modenames{i}(1:end-5);
end

% dot product b/w modes (measure of orthogonality)
dots = modes'*modes;

f = figure;
set(f,'units','normalized','Position',[0.45,0.18,0.55,0.7])
imagesc(dots);
h = colorbar;
h.Label.String = 'Dot Product';
cmap = flip(gray);
colormap(cmap)
caxis([-1 1])
% caxis([min(min(dots)),max(max(dots))])
xticks([1:1+numel(modenames)])
xticklabels(modenames)
yticks([1:1+numel(modenames)])
yticklabels(modenames)
title(titleStr)

ax = gca;
ax.FontSize = 35;

end % dotProductModes