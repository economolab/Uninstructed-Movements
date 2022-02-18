function fn = findDataFnMultiSession(meta,anmix,dateix)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm{anmix}, meta.date{anmix}{dateix}};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');
if iscell(fn)
    fn = fn{1};
end

end % findDataFn