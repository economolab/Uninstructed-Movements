function proj = getProjection(fr, mode)
%fr is time x latents (or cells) x trials
%size of mode should match second dimension of fr

proj = zeros(size(fr, 1), size(fr, 3));     % (time x trials)
for i = 1:size(proj, 2)                     % For each trial...
    proj(:, i) = squeeze(fr(:, :, i))*mode;
end