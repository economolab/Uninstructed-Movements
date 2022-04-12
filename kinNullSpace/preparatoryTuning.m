function tuning_ratio = preparatoryTuning(rez,time,params)

moveix = find(time>=params.move(1) & time<=params.move(2));
prepix = find(time>=params.prep(1) & time<=params.prep(2));

N_null_move = rez.N(moveix,:) * rez.W_null;
N_null_prep = rez.N(prepix,:) * rez.W_null;

N_potent_move = rez.N(moveix,:) * rez.W_potent;
N_potent_prep = rez.N(prepix,:) * rez.W_potent;

gamma = norm(N_null_move,'fro')^2 / ...
        norm(N_potent_move,'fro')^2;


tuning_ratio = (1/gamma) * (norm(N_null_prep,'fro')^2 / norm(N_potent_prep,'fro')^2);


end