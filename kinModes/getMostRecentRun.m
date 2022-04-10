function runix = getMostRecentRun(pth,anm_date)


contents = dir(pth);
% find files that end in .mat
keep = false(size(contents));
for i = 1:numel(contents)
    if contains(contents(i).name,'.mat') && contains(contents(i).name,anm_date)
        keep(i) = true;
    end
end

mats = contents(keep); % just the  .mats

% get file names
fns = {mats.name}';

% sort the filenames
sorted = natsortfiles(fns);

% last file is now most recent. get runix from it
fn = sorted{end};

a = split(fn,{'_','.'});
runix = a{3};

end


