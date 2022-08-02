function meta = loadMAH13_MCStim(meta,datapth)

% % first session, didn't record nidq, not using
% meta(end+1).datapth = datapth;
% meta(end).anm = 'MAH13';
% meta(end).date = '2022-07-28';
% meta(end).datafn = findDataFn(meta(end));
% meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
% meta(end).stim = 'bilateral';
% meta(end).stimLoc = 'MC';

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH13';
meta(end).date = '2022-07-29';
meta(end).datafn = findDataFn(meta(end));
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
meta(end).stim = 'bilateral';
meta(end).stimLoc = 'MC';

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH13';
meta(end).date = '2022-08-01';
meta(end).datafn = findDataFn(meta(end));
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
meta(end).stim = 'bilateral';
meta(end).stimLoc = 'MC';

% meta(end+1).datapth = datapth;
% meta(end).anm = 'MAH13';
% meta(end).date = '2022-08-02';
% meta(end).datafn = findDataFn(meta(end));
% meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
% meta(end).stim = 'bilateral';
% meta(end).stimLoc = 'MC';

end

