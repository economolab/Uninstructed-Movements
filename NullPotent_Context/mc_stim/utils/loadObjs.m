function objs = loadObjs(meta)
objs = cell(1,numel(meta));
for i = 1:numel(meta)
    dat = load(meta(i).datapth);
    objs{i} = dat.obj;
end
end % loadObjs