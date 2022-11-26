function plotHazCDandJaw_multiDell(latent,AvgAll,params,taxis)
for e = 1:3                 % For CDearly, CDlate, and jaw prob...
    for d = 1:4             % For the first 4 delay lengths...
        if e==1             % CDearly
            la = latent.early;  
            dur = 'CDearly';
            s = d;
            colors = {[0 0 1],[1 0 0]};
        elseif e==2         % CDlate
            la = latent.late;
            dur = 'CDlate';
            colors = {[0 0 1],[1 0 0]};
            s = d+4;
        elseif e==3         % Jaw stuff
            la = AvgAll;
            s = d+8;
            colors = {[0 0 0.9],[0.9 0 0]};
        end
        subplot(3,4,s)
        plot(taxis,la.left{d},'Color',colors{2},'LineWidth',2)
        hold on;
        plot(taxis,la.right{d},'Color',colors{1},'LineWidth',2)
        xlim([-1.4 2])
        xline(0,'LineStyle','--','LineWidth',1.3)
        xline(-1.3,'LineStyle',':','LineWidth',1.3)
        xline(params.delay(d),'LineStyle','-.','LineWidth',1.3)
        %legend('Left','Right','Delay onset','Sample onset','GoCue','Location','best')
        len = num2str(params.delay(d));
        if e==1 || e==2 
            subtitle = strcat(dur,';Delay length =',{' '},len);
            ylabel('a.u.')
        elseif e==3
            subtitle = strcat('Jaw prob; Delay length =',{' '},len);
            ylabel('%')
            xlabel('Time since delay onset')
        end
        title(subtitle)
    end
end
end