function meta = loadMAH14_MCStim(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH14';
meta(end).date = '2022-08-15';
meta(end).datafn = findDataFn(meta(end));
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
meta(end).stim = 'bilateral';
meta(end).stimLoc = 'Bi_MC';

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH14';
meta(end).date = '2022-08-16';
meta(end).datafn = findDataFn(meta(end));
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
meta(end).stim = 'bilateral';
meta(end).stimLoc = 'Bi_MC';

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH14';
meta(end).date = '2022-08-17';
meta(end).datafn = findDataFn(meta(end));
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
meta(end).stim = 'bilateral';
meta(end).stimLoc = 'Bi_MC';

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH14';
meta(end).date = '2022-08-18';
meta(end).datafn = findDataFn(meta(end));
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
meta(end).stim = 'bilateral';
meta(end).stimLoc = 'Bi_MC';



end

