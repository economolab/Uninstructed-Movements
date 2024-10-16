% Change cluster quality names
<<<<<<< Updated upstream
obj = deleteGarbageClu(obj);                    % Use function to get rid of all of the clusters that haev 'garbage' label
=======

obj = deleteGarbageClu(obj);                    % Use function to get rid of all of the clusters that have 'garbage' label
>>>>>>> Stashed changes
numprobes = size(obj.clu,2);                    % Number of probes in recording session
for i = 1:numprobes                             % For each probe...
    currprobe = i;
    numclu = size(obj.clu{currprobe},2);        % Get the number of clusters
    for cc = 1:numclu                           % For each cluster...
        qual = obj.clu{currprobe}(cc).quality;  % Get the 'quality' label
        if contains(qual,'good') || contains(qual,'fair') || contains(qual,'poor')  % If the quality label is a good, fair, or poor
            qual = qual(1:4);                   % Take the first four letters of the quality label (get rid of the extraneous spaces)
        elseif contains(qual,'great') || contains(qual,'multi')
            qual = qual(1:5);                   % Take the first five letters of the quality label
        elseif contains(qual,'excellent')
            qual = qual(1:9);                   % Take the first nine letters of the quality label
        end
        obj.clu{currprobe}(cc).quality = qual;  % Save the new quality label to the data object
    end
end
<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes
