tm = obj.time;
f = figure; 
f.Position = [514    54   822   933];
for i = 1:numel(params.trials2use) % trial
    trix = params.trials2use(i);

    clf
    ax1 = subplot(3,1,1);
    hold on;
    imagesc(tm,1:nrDims,squeeze(vidR_Cam1(:,i,:))');
    colorbar;
    caxis([-2 2])
    xlim([-1.3 1])
    ylim([1 200])

    ltm = obj.bp.ev.lickL{trix} - obj.bp.ev.(params.alignEvent)(trix);
    ltm = ltm(ltm>-1 & ltm<1);
    rtm = obj.bp.ev.lickR{trix} - obj.bp.ev.(params.alignEvent)(trix);
    rtm = rtm(rtm>-1 & rtm<1);
    licktm = sort([ltm rtm]);
    hold on;
%     for j = 1:numel(licktm)
%         plot([licktm(j) licktm(j)],[190 200],'k--','LineWidth',2)
%     end
    title('Video SVs')


    hold off;

    ax2 = subplot(3,1,2); sgtitle(['Trial ' num2str(trix)])
    hold on;
    temp = squeeze(absVidR_Cam1(:,i,:));
%     temp = normalizeInRange(temp,[0 10]);
    imagesc(tm,1:nrDims,temp')
    colorbar
    caxis([-4 4])
    xlim([-1.3 1])
    ylim([1 200])

    hold on;
%     for j = 1:numel(licktm)
%         plot([licktm(j) licktm(j)],[190 200],'k--','LineWidth',2)
%     end

    title('ME SVs')
%     pause
    hold off;






    ax3 = subplot(3,1,3); 
    
%     if obj.bp.autowater(i)
%         sgtitle(['Trial ' num2str(trix) ' | Autowater On'])
%     else
%         sgtitle(['Trial ' num2str(trix) ' | Autowater Off'])
%     end
    sgtitle(['Trial ' num2str(trix) ])



    hold on;
    temp = squeeze(N(:,trix,:)); % before trimming trials
%     temp = normalizeInRange(temp,[0 10]);
    imagesc(tm,1:size(temp,2),temp')
    colorbar
%     caxis([0 100])
    xlim([-1.3 1])
    ylim([1 size(temp,2)])

    hold on;
%     for j = 1:numel(licktm)
%         plot([licktm(j) licktm(j)],[size(temp,2)-3 size(temp,2)],'k--','LineWidth',2)
%     end

    title('Neural Data')
    pause
    hold off;


end


%% 1 sv, all trials

tm = obj.time;
f = figure; 
f.Position = [514    54   822   933];
for i = 1:size(vidR,3) % svdims
    clf
    ax1 = subplot(2,1,1);
    hold on;
    imagesc(tm,1:346,squeeze(vidR(:,:,i))');
    colorbar;
    caxis([-4 4])
    xlim([-1 1])
    ylim([1 346])

%     ltm = obj.bp.ev.lickL{i} - obj.bp.ev.(params.alignEvent)(i);
%     ltm = ltm(ltm>-1 & ltm<1);
%     rtm = obj.bp.ev.lickR{i} - obj.bp.ev.(params.alignEvent)(i);
%     rtm = rtm(rtm>-1 & rtm<1);
%     licktm = sort([ltm rtm]);
%     hold on;
%     for j = 1:numel(licktm)
%         plot([licktm(j) licktm(j)],[190 200],'k--','LineWidth',2)
%     end
    title('Video SVs')


    hold off;

    ax2 = subplot(2,1,2); sgtitle(['Trial ' num2str(i)])
    hold on;
    temp = squeeze(absVidR(:,:,i));
%     temp = normalizeInRange(temp,[0 10]);
    imagesc(tm,1:346,temp')
    colorbar
    caxis([-4 4])
    xlim([-1 1])
    ylim([1 346])

%     ltm = obj.bp.ev.lickL{i} - obj.bp.ev.(params.alignEvent)(i);
%     ltm = ltm(ltm>-1 & ltm<1);
%     rtm = obj.bp.ev.lickR{i} - obj.bp.ev.(params.alignEvent)(i);
%     rtm = rtm(rtm>-1 & rtm<1);
%     licktm = sort([ltm rtm]);
%     hold on;
%     for j = 1:numel(licktm)
%         plot([licktm(j) licktm(j)],[190 200],'k--','LineWidth',2)
%     end

    title('ME SVs')

    
    sgtitle(['Dim ' num2str(i)])


%     saveas(f,['Dim ' num2str(i)],'png')
    pause
    hold off;


end




%% first lick histogram

% for i = 1:size(vidR,2) % trial
%     ltm = obj.bp.ev.lickL{i} - obj.bp.ev.(params.alignEvent)(i);
%     ltm = ltm(ltm>0);
%     if ~isempty(ltm)
%         ltm = ltm(1);
%     end
%     rtm = obj.bp.ev.lickR{i} - obj.bp.ev.(params.alignEvent)(i);
%     rtm = rtm(rtm>0);
%     if ~isempty(rtm)
%         rtm = rtm(1);
%     end
%     licktm = sort([ltm rtm]);
%     if ~isempty(licktm)
%         firstLick(i)  = licktm(1);
%     else
%         firstLick(i) = nan;
%     end
% end
% firstLick = firstLick';


