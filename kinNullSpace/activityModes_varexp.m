function varexp = activityModes_varexp(modes,rez)


modeNames = fieldnames(modes.null);


% null space / null activity modes
temp = permute(rez.N_null,[1 3 2]);
temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
[~,~,eigval]= pca(temp);

C = cov(temp);
eigsum = sum(eigval);

for i = 1:numel(modeNames)
    varexp.null_null.(modeNames{i}) = var_proj(modes.null.(modeNames{i}),C,eigsum);
end

% % potent space / null activity modes
% temp = permute(rez.N_potent,[1 3 2]);
% temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
% [~,~,eigval]= pca(temp);
% 
% C = cov(temp);
% eigsum = sum(eigval);
% 
% for i = 1:numel(modeNames)
%     varexp.potent_null.(modeNames{i}) = var_proj(modes.null.(modeNames{i}),C,eigsum);
% end


% % null space / potent activity modes
% temp = permute(rez.N_null,[1 3 2]);
% temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
% [~,~,eigval]= pca(temp);
% 
% C = cov(temp);
% eigsum = sum(eigval);
% 
% for i = 1:numel(modeNames)
%     varexp.null_potent.(modeNames{i}) = var_proj(modes.potent.(modeNames{i}),C,eigsum);
% end

% potent space / potent activity modes
temp = permute(rez.N_potent,[1 3 2]);
temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
[~,~,eigval]= pca(temp);

C = cov(temp);
eigsum = sum(eigval);

for i = 1:numel(modeNames)
    varexp.potent_potent.(modeNames{i}) = var_proj(modes.potent.(modeNames{i}),C,eigsum);
end


end