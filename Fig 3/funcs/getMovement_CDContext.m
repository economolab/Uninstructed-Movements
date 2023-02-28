function rez = getMovement_CDContext(obj,params,cond2use,cond2proj,kin,standardize)

cd_labels = {'context'};
cd_epochs = {'sample'};
cd_times = {[-0.29 -0.1]}; % in seconds, relative to respective epochs

for sessix = 1:numel(obj)

    %-------------------------------------------
    % --setup results struct--
    % ------------------------------------------
    rez(sessix).time = obj(sessix).time;
    if strcmp(standardize,'standardize')
        rez(sessix).movedat = formatKin(kin(sessix).dat_std,params(sessix),rez(sessix).time,cond2use);
    else
        rez(sessix).movedat = formatKin(kin(sessix).dat,params(sessix),rez(sessix).time,cond2use);
    end
    rez(sessix).condition = params(sessix).condition;
    rez(sessix).trialid = params(sessix).trialid;
    rez(sessix).alignEvent = params(sessix).alignEvent;
    rez(sessix).align = mode(obj(sessix).bp.ev.(rez(sessix).alignEvent));
    rez(sessix).ev.sample = obj(sessix).bp.ev.sample;
    rez(sessix).ev.delay = obj(sessix).bp.ev.delay;
    rez(sessix).ev.goCue = obj(sessix).bp.ev.goCue;


    % ------------------------------------------
    % --get coding directions--
    % ------------------------------------------
    rez(sessix).movecd_mode = zeros(size(rez(sessix).movedat,2),numel(cd_labels)); % (movementFeatures,numCDs)
    for ix = 1:numel(cd_labels)
        % find time points to use
        e1 = mode(rez(sessix).ev.(cd_epochs{ix})) + cd_times{ix}(1) - rez(sessix).align;
        e2 = mode(rez(sessix).ev.(cd_epochs{ix})) + cd_times{ix}(2) - rez(sessix).align;
        times.(cd_labels{ix}) = rez(sessix).time>e1 & rez(sessix).time<e2;
        % calculate coding direction
        movedims = [1,2];                       % Which of the cond2use you want to use for finding the mode
        rez(sessix).movecd_mode(:,ix) = calcMovementCD(rez(sessix).movedat,times.(cd_labels{ix}),movedims);
    end

    % ------------------------------------------
    % --orthogonalize coding directions--
    % ------------------------------------------
    rez(sessix).movecd_mode_orth = gschmidt(rez(sessix).movecd_mode);


    % ------------------------------------------
    % --project neural population on CDs--
    % ------------------------------------------
    temp = permute(rez(sessix).movedat(:,:,movedims),[1 3 2]); % (time,cond,movefeatures), permuting to use tensorprod() on next line for the projection
    rez(sessix).movecd_proj = tensorprod(temp,rez(sessix).movecd_mode_orth,3,1); % (time,cond,movecd), cond is in same order as cond2use variable defined at the top of this function


    % ------------------------------------------
    % --variance explained--
    % ------------------------------------------
    movedat = rez(sessix).movedat(:,:,movedims);
    datacov = cov(cat(1,movedat(:,:,1),movedat(:,:,2)));
    datacov(isnan(datacov)) = 0;
    eigsum = sum(eig(datacov));
    for i = 1:numel(cd_labels)
        % whole trial
        rez(sessix).movecd_varexp(i) = var_proj(rez(sessix).movecd_mode_orth(:,i), datacov, eigsum);
        % respective epoch
        epoch_movedat = rez(sessix).movedat(times.(cd_labels{i}),:,movedims);
        epoch_datacov = cov(cat(1,epoch_movedat(:,:,1),epoch_movedat(:,:,2)));
        epoch_datacov(isnan(epoch_datacov)) = 0;
        epoch_eigsum = sum(eig(epoch_datacov));
        rez(sessix).movecd_varexp_epoch(i) = var_proj(rez(sessix).movecd_mode_orth(:,i), epoch_datacov, epoch_eigsum);
    end

    % ------------------------------------------
    % --selectivity--
    % ------------------------------------------
    % coding directions
    rez(sessix).selectivity_squared = squeeze(rez(sessix).movecd_proj(:,1,:) - rez(sessix).movecd_proj(:,2,:)).^2; 
    % sum of coding directions
    rez(sessix).selectivity_squared(:,4) = sum(rez(sessix).selectivity_squared,2); 
    % full neural pop
    temp = rez(sessix).movedat(:,:,movedims);
    temp = (temp(:,:,1) - temp(:,:,2)).^2;
    rez(sessix).selectivity_squared(:,5) = sum(temp,2); % full neural pop

    % ------------------------------------------
    % --selectivity explained--
    % ------------------------------------------
    full = rez(sessix).selectivity_squared(:,5);
    rez(sessix).selexp = zeros(numel(rez(sessix).time),4); % (time,nCDs+1), +1 b/c sum of CDs
    for i = 1:4
        % whole trial
        rez(sessix).selexp(:,i) = 1 - ((full - rez(sessix).selectivity_squared(:,i)) ./ full);
    end


    % set some more rez variables to keep track of
    rez(sessix).cd_times = times;
    rez(sessix).cd_labels = cd_labels;
    rez(sessix).movecd_epochs = cd_epochs;
    rez(sessix).movecd_times_epoch = cd_times; % relative to respective epochs


end

end