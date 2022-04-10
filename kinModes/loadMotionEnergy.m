function me = loadMotionEnergy(meta)

mepth = fullfile(meta.datapth,'DataObjects',meta.anm); 
contents = dir(mepth);
contents = contents(~[contents.isdir]);

fns = {contents.name}';

fn = patternMatchCellArray({contents.name}',{'motionEnergy',meta.date},'all');

if numel(fn) == 1
    fn = fn{1};
else
    disp(['UNABLE TO LOCATE A MOTION ENERGY FILE IN: ' mepth ' . Continuing without motion energy as a feature']);
    me.use = 0; 
    return
end

temp = load(fullfile(mepth,fn));
me = temp.me;
me.use = 1;

end % loadMotionEnergy