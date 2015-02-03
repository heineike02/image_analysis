function norm_val = get_normval(all_tracks_vec, all_times_vec, norm_ind, norm_ph,channel)
%right now this is specific to nf_calcs, but could be for other functions
%that operate on tracks
all_tracks = all_tracks_vec{norm_ind};
all_times = all_times_vec{norm_ind};
tracks = all_tracks.(norm_ph);
timevals = all_times.(norm_ph);
times_ind = 1:length(timevals);
[mean_nf, std_nf] = nf_calcs(times_ind,tracks,channel);

norm_val = mean(mean_nf);