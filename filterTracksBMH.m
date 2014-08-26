function tracks_out = filterTracksBMH(tracks)

nmi_filt = 1;
if nmi_filt == 1;

    %Parameters
    %Const Bright filter
    std_above_mean_nmi = 3;
    %Throw out tracks that have an average nmi greater than nmi_max, defined as
    %std_above_mean_nmi standard deviations above the mean nmi

    if isstruct(tracks(1).nmi)
        
        %Apply filter to each channel
        channels = fields(tracks(1).nmi);
        for ch = 1:length(channels)
            nTracks = length(tracks);
            mean_nmi_vec = zeros(nTracks,1);
            for jj = 1: nTracks
                mean_nmi_vec(jj) = mean([tracks(jj).nmi.(channels{ch})]);
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
            tracks = tracks_out;
        end
            
    else
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
        %reset tracks for next filter
        tracks = tracks_out;

    end

end

%% This filter is for KL_M2MC data from20140416.  Position 3 had one NF point greater than 17 - it might be worth figuring out why that happened and tuning the algorithm
nf_filt = 1;

if nf_filt == 1
    nf_max = 17.0

    %Throw out tracks that have an average nf greater than nf_max,
    
    
    if isstruct(tracks(1).nmi)
        %Apply filter to each channel
        channels = fields(tracks(1).nmi);
        for ch = 1:length(channels)
              
            nTracks = length(tracks);
            mean_nf_vec = zeros(nTracks,1);

            for jj = 1: nTracks
                mean_nf_vec(jj) = mean([tracks(jj).nf.(channels{ch})]);
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
        end
    else
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
    end
end






