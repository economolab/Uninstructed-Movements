function [cd_null, cd_potent, cd_context] = CDContext_AllSpaces_NonStation(obj,meta,rez,popfns,condfns,movefns,MoveNonMove,params,zscored)
for sessix = 1:numel(meta)                                                  % For each session...
    [blockid,nBlocks] = getBlockNum_AltContextTask(sessix,obj);

    for p = 1:length(popfns)                                            % For null, potent, and full population...
        temp1 = []; temp2 = [];
        switch p                                                        % Get the appropriate single trial PSTHs
            case 1
                temppsth = rez(sessix).recon.null;                      % Single-trial PSTHs reconstructed from the null space (time x trials x cells)
            case 2
                temppsth = rez(sessix).recon.potent;                    % Single-trial PSTHs reconstructed from the potent space (time x trials x cells)
            case 3
                %                     temppsth = obj(sessix).trialdat;                        % Single-trial PSTHs from full neural data (time x cells x trials)
                %                     temppsth = permute(temppsth,[1 3 2]);                   % Switch to format (time x trials x cells)
                temppsth = zscored(sessix).trialdat;                        % Single-trial PSTHs from full neural data (time x cells x trials)
        end
 
        for c = 1:length(condfns)                                       % For each condition...
            condtrials = MoveNonMove(sessix).all.(condfns{c});          % Get the trials for this condition
            avgpsth1 = mean(temppsth(:,condtrials,:),2,'omitnan');      % Get the condition-averaged PSTH
            temp1 = cat(2,temp1,avgpsth1);                              % Concatenate both conditions (time x condition x cells)

            for mo = 1:length(movefns)                                                  % For all trials, move trials, and non-move trials...
                cond_MnM_trix = MoveNonMove(sessix).(movefns{mo}).(condfns{c});                 % Get the cond trials
                psth2 = temppsth(:,cond_MnM_trix,:);                                    % Get the single trial PSTHs, separated by condition (time x trials x cells)
                psth2use.(popfns{p}).MnM.(movefns{mo}).(condfns{c}) = permute(psth2,[1 3 2]);    % MnM PSTHs: (time x cells x trials); store separately for each condition and move cond
            end
        end
%         psth2use.(popfns{p}) = temppsth;
    end

    % -- Find coding dimensions from RECONSTRUCTED full neural activity which is reconstructed from the null and potent spaces
    % Use train data to calculate the CD; Project test data onto CD %
    cond2use = [1 2];            % (NUMBERING ACCORDING TO THE CONDITIONS PROJECTED INTO NULL AND POTENT SPACES, i.e. which of the conditions specified in 'cond2proj' above do you want to use?)
    cond2proj = [1 2];           % 2AFC hits/misses, AW hits/misses(corresponding to null/potent psths in rez)
    cond2use_trialdat = [1 2];
    cd_null(sessix) = getCodingDimensions_Context_NonStation(rez(sessix).recon_psth.null,... 
        rez(sessix).recon.null,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,nBlocks,blockid, psth2use.null.MnM);

    cd_potent(sessix) = getCodingDimensions_Context_NonStation(rez(sessix).recon_psth.potent,...
        rez(sessix).recon.potent,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,nBlocks,blockid,psth2use.potent.MnM);

    % Calc CDContext from full neural pop
    cd_context(sessix) = getCodingDimensions_Context_NonStation(obj(sessix).psth,...
        zscored(sessix).trialdat,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,nBlocks,blockid,psth2use.fullpop.MnM);
end