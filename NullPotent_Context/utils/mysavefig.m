function mysavefig(f,pth,fn,forPub)
% f is fig handle
% pth is where to save
% fn is file name    
% if you specify a fourth argument, print figure for publication quality

if ~exist(pth,'dir')
    mkdir(pth)
end

if nargin > 3
    print(f, fullfile(pth,fn), '-depsc', '-r1000')
else
    % savefig(f,fullfile(pth,fn))
    saveas(f,fullfile(pth,fn),'png')
%     saveas(f,fullfile(pth,fn),'svg')
%     saveas(f,fullfile(pth,fn),'epsc')
end




end