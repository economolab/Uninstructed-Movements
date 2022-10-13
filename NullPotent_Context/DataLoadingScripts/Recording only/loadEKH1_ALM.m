function meta = loadEKH1_ALM(meta)
meta(end+1).datapth = fullfile('C:\Users\Jackie\Documents\Grad School\Economo Lab\DataObjects\EKH1');
meta(end).anm = 'EKH1';
meta(end).date = '2021-07-22';
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
meta(end).datapth = fullfile('C:\Users\Jackie\Documents\Grad School\Economo Lab\DataObjects\EKH1');
meta(end).anm = 'EKH1';
meta(end).date = '2021-08-05';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2;

meta(end+1) = meta(end);
meta(end).datapth = fullfile('C:\Users\Jackie\Documents\Grad School\Economo Lab\DataObjects\EKH1');
meta(end).anm = 'EKH1';
meta(end).date = '2021-08-07';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2;

end

function objfn = findDataFn(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');
objfn = fn{1};

end % loadRawDataObj