function objs = loadObjs(meta)
objs = cell(1,numel(meta));
for i = 1:numel(meta)
    dat = load(fullfile(meta(i).datapth, meta(i).datafn));
    objs{i} = dat.obj;
end
end % loadObjs