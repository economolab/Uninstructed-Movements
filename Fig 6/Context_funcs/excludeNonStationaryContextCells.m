function [tempblockpsth, goodcell] = excludeNonStationaryContextCells(sessix, obj, nBlocks,blockid,times)
tempblockpsth = NaN(size(obj(sessix).psth,1),size(obj(sessix).psth,2),nBlocks);     % [time x cells x number of blocks]
for bb = 1:nBlocks                                                                  % For each block number...
    blockix = find(blockid==bb);                                                    % Find the trials in this block
    tempblockpsth(:,:,bb) = mean(obj(sessix).trialdat(:,:,blockix),3,'omitnan');    % Take the avg FR of all cells during this block
end
blockpsth = mean(tempblockpsth(times.start:times.stop,:,:),1,'omitnan');            % Get the avg FR of each cell during the presample period (on each block)
nCells = size(blockpsth,2);
goodcell = NaN(1,nCells);
for cc = 1:nCells                                                                   % For every cell...
    nSwitches = nBlocks-1;                                                          % Number of times the context switches
    switches = NaN(1,nSwitches);                                
    for ss = 1:nSwitches                                                            % For each block switch...
        FRblockA = blockpsth(:,cc,ss);                                              % Get the mean FR for the first block in the pair
        FRblockB = blockpsth(:,cc,ss+1);                                            % Get the mean FR for the subsequent block in the pair
        if FRblockA>FRblockB                                                        % If the FR is higher in the first block
            switches(ss) = 1;                                                       % Note this as 1
        elseif FRblockB>FRblockA                                                    % IF the FR is higher ine the second block
            switches(ss)=0;                                                         % Note this as 0
        end
    end

    
    % A true context-selective cell should alternate between being higher
    % in blockA and blockB depending on which context is blockA
    % A non-stationary cell that is decreasing in FR across the session,
    % for example, will always have a higher FR in block A
    if sum(switches)==(nSwitches/2)                 % If the cell's FR is higher on roughly half of the blocks...
        goodcell(cc) = 1;                               % Count it as a truly context-selective cell   
    elseif sum(switches)==floor(nSwitches/2)
        goodcell(cc) = 1;
    elseif sum(switches)==ceil(nSwitches/2)
        goodcell(cc) = 1;
    elseif sum(switches)==(floor(nSwitches/2)-1)
        goodcell(cc) = 1;
    else                                            % If the cell's selectivity is not consistent across the session...
        goodcell(cc) = 0;                               % Don't count this cell as being truly context-selective
    end
end

    