function [out,kernsd] = mySmooth(x, N)
% operates on first dimension only
% gaussian kernel with window size N, std dev sigma

% returns:
% - out: filtered data
% - kernsd: std dev of gaussian kernel

    Ncol = size(x, 2);
    Nel = size(x, 1);
    
    kern = gausswin(N);
    kernsd = std(1:N);
    
    
    kern(1:floor(numel(kern)/2)) = 0; %causal
    kern = kern./sum(kern);
    
    out = zeros(Nel, Ncol);
    for j = 1:Ncol
        out(:, j) = conv(x(:, j), kern, 'same');
    end

end