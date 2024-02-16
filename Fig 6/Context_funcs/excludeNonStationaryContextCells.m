function [tempblockpsth, goodcell] = excludeNonStationaryContextCells(sessix, obj, nBlocks,blockid,times)
tempblockpsth = NaN(size(obj(sessix).psth,1),size(obj(sessix).psth,2),nBlocks);
for bb = 1:nBlocks
    blockix = find(blockid==bb);
    tempblockpsth(:,:,bb) = mean(obj(sessix).trialdat(:,:,blockix),3,'omitnan');
end
blockpsth = mean(tempblockpsth(times.start:times.stop,:,:),1,'omitnan');
nCells = size(blockpsth,2);
goodcell = NaN(1,nCells);
for cc = 1:nCells
    nSwitches = nBlocks-1;
    switches = NaN(1,nSwitches);
    for ss = 1:nSwitches
        FRblockA = blockpsth(:,cc,ss);
        FRblockB = blockpsth(:,cc,ss+1);
        if FRblockA>FRblockB
            switches(ss) = 1;
        elseif FRblockB>FRblockA
            switches(ss)=0;
        end
    end
    %         if sum(switches)==nSwitches
    %             goodcell(cc) = 0;
    %         elseif sum(switches)==0
    %             goodcell(cc) = 0;
    %         elseif sum(switches)==1||sum(switches)==nSwitches-1
    %             goodcell(cc) = 0;
    %         elseif sum(switches)==(nSwitches/2)||sum(switches)==(nSwitches/2)-1||sum(switches)==(nSwitches/2)+1
    %             goodcell(cc) = 1;
    %         end

    if sum(switches)==(nSwitches/2)
        goodcell(cc) = 1;
    elseif sum(switches)==floor(nSwitches/2)
        goodcell(cc) = 1;
    elseif sum(switches)==ceil(nSwitches/2)
        goodcell(cc) = 1;
    elseif sum(switches)==(floor(nSwitches/2)-1)
        goodcell(cc) = 1;
    else
        goodcell(cc) = 0;
    end
end

    