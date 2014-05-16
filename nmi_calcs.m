function [mean_nmi, std_nmi] = nmi_calcs(times,tracks)
nTimes = length(times);
nTracks = length(tracks);

mean_nmi = zeros(nTimes,1);
std_nmi = zeros(nTimes,1);
for kk = 1:nTimes
    time = times(kk);
    nmi_vec = zeros(nTracks,1);
    for jj = 1:nTracks
        track_times = [tracks(jj).times];
        track_nmis = [tracks(jj).nmi];
        time_match = (track_times==time);
        if sum(time_match) ~= 0
            nmi_vec(jj) = track_nmis(time_match);
        end
    end
    nmi_vec = nmi_vec(nmi_vec>0);
    mean_nmi(kk) = mean(nmi_vec);
    std_nmi(kk) = std(nmi_vec);
end

end
