function filtered_tracks_cell = filter_singlecell_Analysis(tracks, field, Xper)

    % convert structure to cell
    tracks_cell(:,:) = struct2cell(tracks);
    % isolate the single cell traces
    desired_track_cell = tracks_cell(field,:);
    % access the .RFP (.channel) info
    desired_track_cell = cell2mat(desired_track_cell);
    % access the .RFP
    temp_track_var(:,:) = struct2cell(desired_track_cell);
    desired_track_cell = temp_track_var;
    % obtain the lengths of each trace
    track_lengths = cellfun(@length, desired_track_cell, 'UniformOutput', false);
    track_lengths = cell2mat(track_lengths);
    % only save the traces that are >= 50% of the maximum length
    XPercent = max(track_lengths)*Xper;
    indc_of_tracks = find(track_lengths >= XPercent);
    filtered_tracks_cell = tracks_cell(:,indc_of_tracks);