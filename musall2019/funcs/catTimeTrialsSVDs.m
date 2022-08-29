function catted = catTimeTrialsSVDs(dat)

temp = cell(size(dat));

for i = 1:numel(temp)
    temp{i} = reshape(dat{i},size(dat{i},1)*size(dat{i},2),size(dat{i},3));
end
catted = cat(1,temp{:}); % concatenate conditions

end