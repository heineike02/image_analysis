function [mean_nf, std_nf] = nf_calcs(timevals,tracks,channel)

nTimes = length(timevals);
nTracks = length(tracks); 

mean_nf = zeros(nTimes,1);
std_nf = zeros(nTimes,1);
for kk = 1:nTimes
    timeval = timevals(kk);
    nf_vec = zeros(nTracks,1);
    for jj = 1:nTracks
        track_times = [tracks(jj).times];
        track_nfs = [tracks(jj).nf.(channel)];
        time_match = (track_times==timeval);
        if sum(time_match) ~= 0
            nf_vec(jj) = track_nfs(time_match);
        end
        nf_vec = nf_vec(nf_vec>0);
        
        % median and interquartile range
        mean_nf(kk) = median(nf_vec);
        std_nf(kk) = iqr(nf_vec);
        % mean and standard deviation
        %mean_nf(kk) = mean(nf_vec);
        %std_nf(kk) = std(nf_vec);
    end

end
end
