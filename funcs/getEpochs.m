function [prepix, moveix] = getEpochs(time, prepEpoch, moveEpoch)

prepix = getEpochIx(time,prepEpoch);
moveix = getEpochIx(time,moveEpoch);

end % getEpochs

%% 
function ix = getEpochIx(time,epoch)
    ix = zeros(size(epoch));
    for i = 1:numel(epoch)
        [~,ix(i)] = min(abs(time-epoch(i)));
    end
    ix = ix(1):ix(2);
end % getEpochIx
