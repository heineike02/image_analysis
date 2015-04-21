function [mean_nmi, std_nmi] = nmi_calcs(timevals,tracks,channel)
nTimes = length(timevals);
nTracks = length(tracks);

mean_nmi = zeros(nTimes,1);
std_nmi = zeros(nTimes,1);
for kk = 1:nTimes
    timeval = timevals(kk);
    nmi_vec = zeros(nTracks,1);
    for jj = 1:nTracks
        track_times = [tracks(jj).times];
        if isempty(channel)
            track_nmis = [tracks(jj).nmi];
        else
            track_nmis = [tracks(jj).nmi.(channel)];
        end
        time_match = (track_times==timeval);
        if sum(time_match) ~= 0
            nmi_vec(jj) = track_nmis(time_match);
        end
    end
    nmi_vec = nmi_vec(nmi_vec>0);
    mean_nmi(kk) = mean(nmi_vec);
    std_nmi(kk) = std(nmi_vec);
end

end
