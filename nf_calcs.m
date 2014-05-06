function [mean_nf, std_nf] = nf_calcs(times,tracks)
nTimes = length(times);
nTracks = length(tracks);

mean_nf = zeros(nTimes,1);
std_nf = zeros(nTimes,1);
for kk = 1:nTimes
    time = times(kk);
    nf_vec = zeros(nTracks,1);
    for jj = 1:nTracks
        track_times = [tracks(jj).times];
        track_nfs = [tracks(jj).nf];
        time_match = (track_times==time);
        if sum(time_match) ~= 0
            nf_vec(jj) = track_nfs(time_match);
        end
    end
    nf_vec = nf_vec(nf_vec>0);
    mean_nf(kk) = mean(nf_vec);
    std_nf(kk) = std(nf_vec);
end

end
