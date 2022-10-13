function meta = loadJEB14_M1TJVideo(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'JEB14';
meta(end).date = '2022-08-22';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2; 
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);


meta(end+1).datapth = datapth;
meta(end).anm = 'JEB14';
meta(end).date = '2022-08-23';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2; 
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'JEB14';
meta(end).date = '2022-08-24';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2; 
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);



end

