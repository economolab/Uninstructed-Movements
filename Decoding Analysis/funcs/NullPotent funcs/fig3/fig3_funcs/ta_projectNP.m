function rez = ta_projectNP(input_data,rez,cond2use,params)

dims = size(input_data);
toproj = reshape(input_data,dims(1)*dims(2),dims(3));

nullproj = toproj * rez.Qnull;
rez.N_null = reshape(nullproj,dims(1),dims(2),rez.dPrep);
potentproj = toproj * rez.Qpotent;
rez.N_potent = reshape(potentproj,dims(1),dims(2),rez.dMove);

nullproj = reshape(nullproj,dims(1),dims(2),rez.dPrep);
potentproj = reshape(potentproj,dims(1),dims(2),rez.dMove);

rez.N_potent_psth = zeros(dims(1),rez.dMove,numel(cond2use));
rez.N_null_psth   = zeros(dims(1),rez.dPrep,numel(cond2use));

for i = 1:numel(cond2use)
    trix = params.trialid{cond2use(i)};

    rez.N_potent_psth(:,:,i) = squeeze(mean(potentproj(:,trix,:),2));
    rez.N_null_psth(:,:,i) = squeeze(mean(nullproj(:,trix,:),2));
end


end