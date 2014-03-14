function tracks_out = filterTracksBMH(tracks)

%Parameters
%Const Bright filter
std_above_mean_nmi = 3;
%Throw out tracks that have an average nmi greater than nmi_max, defined as
%std_above_mean_nmi standard deviations above the mean nmi

nTracks = length(tracks);
mean_nmi_vec = zeros(nTracks,1);

for jj = 1: nTracks
    mean_nmi_vec(jj) = mean([tracks(jj).nmi]);
end

mean_nmi_std = std(mean_nmi_vec);
mean_nmi_mean = mean(mean_nmi_vec);

nmi_max = mean_nmi_mean+std_above_mean_nmi*mean_nmi_std;

nmi_filt = [];
for jj = 1:nTracks
    if mean_nmi_vec(jj)>nmi_max
        nmi_filt = [nmi_filt,jj];
        ['threw out track ', int2str(jj), ' because mean nmi of ',num2str(mean_nmi_vec(jj)),' was greater than nmi_max of ', num2str(nmi_max)]
    end
end




tracks_out = tracks;
tracks_out(nmi_filt) = [];




