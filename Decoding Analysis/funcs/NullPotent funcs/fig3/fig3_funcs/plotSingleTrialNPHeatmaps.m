function plotSingleTrialNPHeatmaps(rez,params,me,ndims,cond2use,meta)

cmap = parula;

sm = 31;


% - null
for sessix = 1:numel(rez)
    f = figure;
    f.Position = [59   139   743   824];
    
    % trials
    trix = [];
    for condix = 1:numel(cond2use)
        trix = [trix ; params(sessix).trialid{cond2use(condix)}];
        condend(condix) = numel(trix);
    end

    % sort dimensions by variance explained (most to least)
    [~,ix] = sort(rez(sessix).ve.null,'descend');

    
    temp = rez(sessix).N_null(:,trix,ix);
    move = me(sessix).move(:,trix);
    medat = me(sessix).data(:,trix);
    medat(isnan(medat)) =0;
    
    for dimix = 1:ndims
        ax = nexttile;
%         toplot = temp(:,:,dimix);
        toplot = temp(:,:,dimix) .* move;
        imagesc(mySmooth(toplot,sm)'); colormap(cmap)
        c = colorbar;
        c.Limits = c.Limits ./ 3;
        for i = 1:numel(condend)
            yline(condend(i),'k--')
        end

%         ax = nexttile;
%         toplot2 = toplot;
%         toplot2(move) = -max(max(toplot));
%         imagesc(mySmooth(toplot2,sm)'); colormap(cmap)
%         c = colorbar;
%         c.Limits = c.Limits ./ 2;
% 
%         ax = nexttile;
%         %   toplot3 = toplot .* ~move;
%         toplot3 = toplot;
%         toplot3(~move) = -max(max(toplot));
%         imagesc(mySmooth(toplot3,sm)'); colormap(cmap)
%         c = colorbar;
%         c.Limits = c.Limits ./ 2;
    end
    sgtitle(['Null | ' meta(sessix).anm ' - ' meta(sessix).date])
end

% - potent
for sessix = 1:numel(rez)
    f = figure;
    f.Position = [992   145   743   824];

    % trials
    trix = [];
    for condix = 1:numel(cond2use)
        trix = [trix ; params(sessix).trialid{cond2use(condix)}];
        condend(condix) = numel(trix);
    end

    % sort dimensions by variance explained (most to least)
    [~,ix] = sort(rez(sessix).ve.potent,'descend');

    
    temp = rez(sessix).N_potent(:,trix,ix);
    move = me(sessix).move(:,trix);
    medat = me(sessix).data(:,trix);
    medat(isnan(medat)) =0;
    
    for dimix = 1:ndims
        ax = nexttile;
%         toplot = temp(:,:,dimix);
        toplot = temp(:,:,dimix) .* move;
        imagesc(mySmooth(toplot,sm)'); colormap(cmap)
        c = colorbar;
        c.Limits = c.Limits ./ 3;
        for i = 1:numel(condend)
            yline(condend(i),'k--')
        end

%         ax = nexttile;
%         toplot2 = toplot;
%         toplot2(move) = -max(max(toplot));
%         imagesc(mySmooth(toplot2,sm)'); colormap(cmap)
%         c = colorbar;
%         c.Limits = c.Limits ./ 2;
% 
%         ax = nexttile;
%         %   toplot3 = toplot .* ~move;
%         toplot3 = toplot;
%         toplot3(~move) = -max(max(toplot));
%         imagesc(mySmooth(toplot3,sm)'); colormap(cmap)
%         c = colorbar;
%         c.Limits = c.Limits ./ 2;
    end
    sgtitle(['Potent | ' meta(sessix).anm ' - ' meta(sessix).date])
end


end
