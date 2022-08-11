function rez = kinematicNullAndPotentSpace(obj,input_data,dfparams,params,me)

clear newrez
newrez = estimateW(obj,input_data,dfparams,obj.time,me); % N,V are zscored neural activity and feature matrix



% NULL AND POTENT SPACE OF W

[newrez.W_null,newrez.W_potent,newrez.N_null,newrez.N_potent,newrez.dPrep,newrez.dMove] = ...
    getNullPotentSpaces_SVD(newrez.W',newrez,input_data); % transposing W b/c we found W by solving V'=N'*W, rather than V = WN

% correlation between V and V_hat = rez.N * rez.W
% rez.corrcoef = calcCorrCoef(rez.V,rez.N,rez.W);

% variance explained by null and potent space
verez = calVarExp_kin(newrez);

rezfns = fieldnames(verez);
for j = 1:numel(rezfns)
    rez.(rezfns{j}) = verez.(rezfns{j});
end


end