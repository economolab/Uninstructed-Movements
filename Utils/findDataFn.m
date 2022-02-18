function objfn = findDataFn(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');
objfn = fn{1};

end % loadRawDataObj

