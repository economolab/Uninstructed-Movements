function rez = getCDs_wRamping(input_data,trialdat_zscored,obj,params,cond2use,cond2use_trialdat,cond2proj,rampcond)
%-----------------------------------------------------------------------------------------------------------------------------
% Slightly different than Nuo Li's version ('Thalamus-driven functional populations in frontal cortex support decision-making')
% Nuo's version: (delay - pre-sample) where delay = population response during the last 500 ms of the delay period and pre-sample = pop response
% during the last 500 ms of pre-sample period 
% This version: (late delay - late sample period) where delay = pop response during last 500 ms of delay period and late sample = last 400 ms of sample period
% Trials used: All R and L correct trials
%-------------------------------------------------------------------------------------------------------------------------

cd_labels = {'early','late','go','ramping'};
cd_epochs = {'delay','goCue','goCue'};
cd_times = {[-0.42 -0.02], [-0.42 -0.02], [0.02 0.42]}; % in seconds, relative to respective epochs
ramp_epochs = {'delay','goCue'};
%ramp_times = {[-0.4 -0.02], [-0.5 -0.02]}; % in seconds, relative to respective epochs
ramp_times = {[-0.6 -0.02], [-0.5 -0.02]}; % in seconds, relative to respective epochs




%-------------------------------------------
% --setup results struct--
% ------------------------------------------
rez.time = obj.time;
rez.psth = input_data;
rez.condition = params.condition;
rez.trialid = params.trialid;
rez.alignEvent = params.alignEvent;
rez.align = mode(obj.bp.ev.(rez.alignEvent));
rez.ev.sample = obj.bp.ev.sample;
rez.ev.delay = obj.bp.ev.delay;
rez.ev.goCue = obj.bp.ev.goCue;


% ------------------------------------------
% --get coding directions--
% ------------------------------------------
rez.cd_mode = zeros(size(rez.psth,2),numel(cd_labels)); % (neurons,numCDs)
for ix = 1:numel(cd_labels)
    if strcmp(cd_labels{ix},'ramping')
        % find time points to use
        e1_start = mode(rez.ev.(ramp_epochs{1})) + ramp_times{1}(1) - rez.align;
        e1_stop = mode(rez.ev.(ramp_epochs{1})) + ramp_times{1}(2) - rez.align;
        times.ramp_lateSamp = rez.time>e1_start & rez.time<e1_stop;

        e2_start = mode(rez.ev.(ramp_epochs{2})) + ramp_times{2}(1) - rez.align;
        e2_stop = mode(rez.ev.(ramp_epochs{2})) + ramp_times{2}(2) - rez.align;
        times.ramp_lateDel = rez.time>e2_start & rez.time<e2_stop;
        % calculate ramping mode
        rez.cd_mode(:,ix) = calcRampingMode(rez.psth,times,rampcond);
    else
        % find time points to use
        e1 = mode(rez.ev.(cd_epochs{ix})) + cd_times{ix}(1) - rez.align;
        e2 = mode(rez.ev.(cd_epochs{ix})) + cd_times{ix}(2) - rez.align;
        times.(cd_labels{ix}) = rez.time>e1 & rez.time<e2;
        % calculate coding direction
        rez.cd_mode(:,ix) = calcCD(rez.psth,times.(cd_labels{ix}),cond2use);
    end
end


% ------------------------------------------
% --orthogonalize coding directions--
% ------------------------------------------
rez.cd_mode_orth = gschmidt(rez.cd_mode);


% ------------------------------------------
% --project neural population on CDs--
% ------------------------------------------
temp = permute(rez.psth(:,:,cond2proj),[1 3 2]); % (time,cond,neurons), permuting to use tensorprod() on next line for the projection
rez.cd_proj = tensorprod(temp,rez.cd_mode_orth,3,1); % (time,cond,cd), cond is in same order as con2use variable defined at the top of this function

% % single trial n/p projs
% proj = reshape(np_trialdat,size(np_trialdat,1)*size(np_trialdat,2),size(np_trialdat,3)) * rez.cd_mode_orth;
% rez.trialdat = reshape(proj,size(np_trialdat,1),size(np_trialdat,2),size(rez.cd_mode_orth,2));


% ------------------------------------------
% --variance explained--
% ------------------------------------------
% psth = rez.psth(:,:,cond2use);
% datacov = cov(cat(1,psth(:,:,1),psth(:,:,2)));
% datacov(isnan(datacov)) = 0;
% eigsum = sum(eig(datacov));
% 
% for i = 1:numel(cd_labels)
%     % whole trial
%     rez.cd_varexp(i) = var_proj(rez.cd_mode_orth(:,i), datacov, eigsum);
%     % respective epoch
%     epoch_psth = rez.psth(times.(cd_labels{i}),:,cond2use);
%     epoch_datacov = cov(cat(1,epoch_psth(:,:,1),epoch_psth(:,:,2)));
%     epoch_datacov(isnan(epoch_datacov)) = 0;
%     epoch_eigsum = sum(eig(epoch_datacov));
%     rez.cd_varexp_epoch(i) = var_proj(rez.cd_mode_orth(:,i), epoch_datacov, epoch_eigsum);
% end

% ------------------------------------------
% --selectivity--
% ------------------------------------------
% coding directions
rez.selectivity_squared = squeeze(rez.cd_proj(:,1,:) - rez.cd_proj(:,2,:)).^2;
% sum of coding directions
rez.selectivity_squared(:,4) = sum(rez.selectivity_squared,2);
% full neural pop
% temp = rez.psth(:,:,cond2use);
% temp = (temp(:,:,1) - temp(:,:,2)).^2;
% temp = (obj.psth(:,:,2) - obj.psth(:,:,3)).^2;
temp = squeeze(mean(trialdat_zscored(:,params.trialid{cond2use_trialdat(1)},:),2)) - squeeze(mean(trialdat_zscored(:,params.trialid{cond2use_trialdat(2)},:),2));
rez.selectivity_squared(:,5) = sum(temp.^2,2); % full neural pop


% ------------------------------------------
% --selectivity explained--
% ------------------------------------------
full = rez.selectivity_squared(:,5);
rez.selexp = zeros(numel(rez.time),4); % (time,nCDs+1), +1 b/c sum of CDs
for i = 1:4
    % whole trial
    rez.selexp(:,i) = rez.selectivity_squared(:,i) ./ full;
end


% set some more rez variables to keep track of
rez.cd_times = times;
rez.cd_labels = cd_labels;
rez.cd_epochs = cd_epochs;
rez.cd_times_epoch = cd_times; % relative to respective epochs

% % quick plot of cds
% for i = 1:3
% figure;
% hold on;
% plot(obj.time,rez.cd_proj(:,1,i),'b')
% plot(obj.time,rez.cd_proj(:,2,i),'r')
% end


end