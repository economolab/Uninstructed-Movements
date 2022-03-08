function analyzeJaw(obj)


%% params
% 
% figure(1)
% set(gcf,'units','normalized','Position',[0.2 0.1 0.7 0.8]); 


dt = 1/400;
sm = 51; % smoothing window

%%

figure; hold on
ct = 0;
for trix = 1:obj.bp.Ntrials
    
    if obj.bp.R(trix)&obj.bp.hit(trix)&~obj.bp.early(trix)&~obj.bp.autowater(trix)
        clr = 'b';
    elseif obj.bp.L(trix)&obj.bp.hit(trix)&~obj.bp.early(trix)&~obj.bp.autowater(trix)
        clr = 'r';
    else 
        continue
    end
    
    tm = obj.traj{1}(trix).frameTimes - 0.5;
    x = mySmooth(obj.traj{1}(trix).ts(:,1,2),sm);
    y = mySmooth(obj.traj{2}(trix).ts(:,2,8),sm);
    z = mySmooth(obj.traj{1}(trix).ts(:,2,2),sm);
    
    tm = tm - obj.bp.ev.jawOnset(trix);
    
    plot(tm(100:end),y(100:end),clr)
    
    ct = ct + 1;
    
end
xlim([-1,0.5])
hold off

end % analyzeJaw




















