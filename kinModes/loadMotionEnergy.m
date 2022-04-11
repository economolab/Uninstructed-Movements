function me = loadMotionEnergy(obj,meta,params,trialnums)

% find and load motion energy .mat file
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
me.use = 1; % if a file was found and successfully loaded, set me.use to 1

% trim number of trials
me.data = me.data(trialnums);

% trim trial length (me.data contains motion energy for each time point in
% trial at 400 Hz). Want to align to params.alignEvent and want to put it
% in same dt as neural data
alignTimes = obj.bp.ev.(params.alignEvent)(trialnums);
me.newdata = zeros(numel(obj.time),numel(trialnums));
for i = 1:numel(trialnums)
    me.newdata(:,i) = interp1(obj.traj{1}(trialnums(i)).frameTimes-0.5-alignTimes(i),me.data{i},obj.time);
end

% replace me.data with me.newdata
me.data = me.newdata;
me = rmfield(me,'newdata');

end % loadMotionEnergy













