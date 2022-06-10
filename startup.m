
dirs = {'DataLoadingScripts',
        'fig1',
        'fig2',
        'fig3',
        'fig4_5',
        'funcs',
        'kinNullSpace',
        'optimization',
        'utils'};
    
for i = 1:numel(dirs)
    addpath(genpath(fullfile(pwd,dirs{i})));
end

