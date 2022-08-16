function psth = standardizePSTH(obj)
% standardize by subtracting baseline firing rate, dividing by baseline std dev

psth = zeros(size(obj.psth));
for i = 1:size(obj.psth,3)
    
    temp = obj.psth(:,:,i);
    mu = obj.presampleFR(:,i);
%     mu = mean(temp)';
    sd = obj.presampleSigma(:,i);
%     sd = std(temp)';
    
%     psth(:,:,i) = (temp - mu') ./ sd';
    psth(:,:,i) = temp;
    
end

end