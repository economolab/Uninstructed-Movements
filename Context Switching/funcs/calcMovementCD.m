function cd = calcMovementCD(movedat,times,movedims)
% calculate a coding direction given:
% movedat - (time,movement features,conditions)
% times - indices or logical array of size(psth,1)=time
%           specifies time points to use when calculating CD
% movedims - which of size(movedat,3)=conditions to use when calculating CD

    tempdat = movedat(:,:,movedims);
    mu = squeeze(mean(tempdat(times,:,:),1));

    sd = squeeze(std(tempdat(times,:,:),[],1));

    cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
%     cd = cd ./ norm(cd);
    cd = cd./sum(abs(cd)); % (ncells,1)
end