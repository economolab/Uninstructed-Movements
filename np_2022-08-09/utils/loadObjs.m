function obj = loadObjs(meta)
for i = 1:numel(meta)
    dat = load(meta(i).datapth);
    if isfield(dat.obj,'meta')
        dat.obj = rmfield(dat.obj,'meta');
    end
    if isfield(dat.obj,'ex')
        dat.obj = rmfield(dat.obj,'ex');
    end
    obj(i) = dat.obj;
    clear dat
end
end % loadObjs