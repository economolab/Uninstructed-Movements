function meta = loadJEB4_ALM(meta)
meta(end+1).datapth = fullfile('C:\Users\Jackie\Documents\Grad School\Economo Lab\DataObjects\JEB4');
meta(end).anm = 'JEB4';
meta(end).date = '2021-03-14';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;

% analysis meta data
meta(end).tmin = -2.5; % (s) relative to params.alignEvent
meta(end).tmax = 3;  % (s) relative to params.alignEvent
meta(end).dt = 0.005;
meta(end).smooth = 15; % smooth psth
% clusters (these qualities are included)
meta(end).quality = {'Fair','Good','Great','Excellent','single'}; 

taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);

% use most of the same fields across sessions
meta(end+1) = meta(end);
meta(end).datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\DataObjects\JEB4';
meta(end).anm = 'JEB4';
meta(end).date = '2021-03-14';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2;

meta(end+1) = meta(end); 
meta(end).datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\DataObjects\JEB4';
meta(end).anm = 'JEB4';
meta(end).date = '2021-03-15';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2;

meta(end+1) = meta(end); 
meta(end).datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\DataObjects\JEB4';
meta(end).anm = 'JEB4';
meta(end).date = '2021-03-16';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;


end

function objfn = findDataFn(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');
objfn = fn{1};

end % loadRawDataObj