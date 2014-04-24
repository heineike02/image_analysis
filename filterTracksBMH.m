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

%% This filter is for KL_M2MC data from20140416.  Position 3 had one NF point greater than 17 - it might be worth figuring out why that happened and tuning the algorithm
nf_max = 17.0

%Throw out tracks that have an average nf greater than nf_max,

nTracks = length(tracks);
mean_nf_vec = zeros(nTracks,1);

for jj = 1: nTracks
    mean_nf_vec(jj) = mean([tracks(jj).nf]);
end

nf_filt = [];
for jj = 1:nTracks
    if mean_nf_vec(jj)>nf_max
        nf_filt = [nf_filt,jj];
        ['threw out track ', int2str(jj), ' because mean nf of ',num2str(mean_nf_vec(jj)),' was greater than nf_max of ', num2str(nf_max)]
    end
end

tracks_out = tracks;
tracks_out(nf_filt) = [];





