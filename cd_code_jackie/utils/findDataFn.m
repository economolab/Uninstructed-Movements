function fn = findDataFn(meta)
objpth = fullfile(meta.datapth,'DataObjects',meta.anm); 
contents = dir(objpth);
contents = contents(~[contents.isdir]);

fns = {contents.name}';

fn = patternMatchCellArray({contents.name}',{'data_structure',meta.date},'all');

if numel(fn) == 1
    fn = fn{1};
else
    error(['Found multiple data_structures containing strings: ' meta.anm ' and/or ' meta.date ' in: ' objpth]);
end

end