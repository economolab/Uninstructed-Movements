function dat = HDF5ToStruct(fpath,fname)

% view contents of hdf5 file
% h5disp(fullfile(fpath,fname));

% get file hierarchy
info = h5info(fullfile(fpath,fname));
%       Filename: '/Users/Munib/Documents/Economo-Lab/code/lfads-jack/model_runs__valid_posterior_sample_and_average'
%           Name: '/'
%         Groups: []
%       Datasets: [12×1 struct]
%      Datatypes: []
%          Links: []
%     Attributes: []

% this only works if the data is organized into datasets (no groups)
for i = 1:numel(info.Datasets)
    dsetNames{i} = info.Datasets(i).Name;
    dat.(dsetNames{i}) = h5read(fullfile(fpath,fname),['/' dsetNames{i}]);
end

end