function psth = normalizePSTH(obj)
% normalize by subtracting baseline firing rate, dividing by baseline std dev

psth = zeros(size(obj.psth));
for i = 1:size(obj.psth,3)
    
    temp = obj.psth(:,:,i);
    mu = obj.presampleFR(:,i);
    sd = obj.presampleSigma(:,i);
    
    psth(:,:,i) = (temp - mu') ./ sd';
    
end

end