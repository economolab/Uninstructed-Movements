function meta = assignDataPath(meta)

for i = 1:numel(meta)
    if ispc
        meta(i).datapth = 'C:\Code\uninstructedMovements-Munib/data';
    else
        meta(i).datapth = '/Users/Munib/Documents/Economo-Lab/data/';
    end
end

end