function me = stitchAndPurgeMoveBouts(me,params,stitch_dist,purge_dist)
% written by jack vincent
% adapted/ported to matlab by munib hasnain

fs = 1 ./ params(1).dt;
stitch_dist = round(stitch_dist*fs);
purge_dist = round(purge_dist*fs);

for sessix = 1:numel(me)
    for trix = 1:size(me(sessix).data,2)

        medat = mySmooth(me(sessix).data(:,trix),11); % for plotting and sanity checks
        medat = fillmissing(medat,"nearest");
        move = me(sessix).move(:,trix);
        move = fillmissing(move,"nearest");

        binary_mask = move;

        % add zero padding at end of mask
        binary_mask = [binary_mask ; zeros(stitch_dist,1)];

        % iterate up to where padding was added
        for i = 1:numel(move) % numel(move)
            % check for falling edge
            if move(i) == 1 && binary_mask(i+1) == 0 % move(i) == 1 && binary_mask(i+1)
                % if a 1 appears anywhere after the falling edge within stitching
                % distance, stitch the 1s together (stitch to the farthest away 1)
                if any(binary_mask(i+1:i+stitch_dist))
                    % find index of last 1 in binary_mask(i+1:i+stitch_dist)
                    last_true = find(binary_mask(i+1:i+stitch_dist),1,'last');
                    binary_mask(i+1:i+1+last_true) = 1;
                end
            end
        end

         % do stitching again
        for i = 1:numel(move) % numel(move)
            % check for falling edge
            if binary_mask(i) == 1 && binary_mask(i+1) == 0 % move(i) == 1 && binary_mask(i+1)
                % if a 1 appears anywhere after the falling edge within stitching
                % distance, stitch the 1s together (stitch to the farthest away 1)
                if any(binary_mask(i+1:i+stitch_dist))
                    % find index of last 1 in binary_mask(i+1:i+stitch_dist)
                    last_true = find(binary_mask(i+1:i+stitch_dist),1,'last');
                    binary_mask(i+1:i+1+last_true) = 1;
                end
            end
        end

        % reset the zero padding just in case
        binary_mask = [binary_mask ; zeros(stitch_dist,1)];

        % will hold locations of rising and falling edges
        rising_idx = [];
        falling_idx = [];

        % in case signal starts high, because that wouldn't be picked up as a rising
        % edge, but it should be
        if binary_mask(1) == 1
            rising_idx(1) = 1;
        end

        % iterate up to where the padding was added
        for i = 1:numel(move) % numel(move)
            % detect rising edges
            if binary_mask(i) == 0 && binary_mask(i+1) == 1 % binary_mask(i) == 1 && binary_mask(i+1) == 1
                rising_idx = [rising_idx ; i+1];
                % detect falling edges
            elseif binary_mask(i) == 1 && binary_mask(i+1) == 0 % binary_mask(i) == 1 && binary_mask(i+1) == 0
                falling_idx = [falling_idx ; i];
            end
        end

        if numel(rising_idx) == (numel(falling_idx) + 1)
            falling_idx = [falling_idx ; numel(move)]; % numel(move)
        end

        % iterate through putative epochs
        for i = 1:numel(rising_idx)
            % the rising idx corresponds to the first 1, the falling idx to the last,
            % so you need to add 1 to their difference to get the total epoch length
            epoch_len = falling_idx(i) - rising_idx(i) + 1;
            % putative epoch length must exceed purge distance, or it gets purged
            if epoch_len < purge_dist
                binary_mask(rising_idx(i):falling_idx(i)+1) = 0;
            end
        end

        % remove zero padding at end
        binary_mask = binary_mask(1:numel(move));

        me(sessix).move(:,trix) = binary_mask;

    end
end


% 
%     clf
%     subplot(2,1,1)
%     hold on;
%     plot(medat)
%     plot(move*120)
% 
%     subplot(2,1,2)
%     hold on;
%     plot(medat)
%     plot(binary_mask*120)
% 
%     pause


end

























