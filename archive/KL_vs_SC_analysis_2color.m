function [all_tracks, all_times] = KL_vs_SC_analysis_2color(ipdir,storeim,fname_conv,op_amp,std_threshSP,species,imdir, maxdisp_1x,pos_fnames,channels,channel_to_image,time_calc,imbg)
% Code for quantifying KL or SC nuclear localization using 1x or 1.5x
% optical zoom.  
% Change jacobs image code output to sort in image order. 
% Too many tracks - something wierd going on
% Had to use ".tiff" for several images
% Note: to make training image use ginput function

'Try to use KL_vs_SC_analysis.m instead of KL_vs_SC_analysis_2color'

%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number, then that is interpreted as a dt between images.  If
%the input is a filename, then that is interpreted as metadata that gives
%exact time values.  

%profile on
%ipdir: directory for image processing
%Windows: 

%Mac
%ex: imdir = '/Users/ucsf/Google Drive/UCSF/ElSamad_Lab/PKA/WetLab/20130716/Gluc_pre_ind/Pos1/'

%fname_conv
% Filename convention
% 'JSO' or 'Micromanager'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files

% op_amp
% optical amplification 
% 1x or 1.5x

if strcmp(op_amp,'1x')
    %load image of KLactis cell - this was made in BMH_20130916_analysis_klac.m
    load([ipdir, 'circKL_1x.mat'])
    circSP.KL = circ;
    %load image of S.Cerevisiae cell (made by Jacob)
    load([ipdir,'circSC_1x.mat'])
    circSP.SC = circ;
    
    %Size of image of individual cell to analyze
    sizSP.KL = [15,15];
    sizSP.SC = [17,17];
    %radius of cell for nuclear enrichment calculation
    radSP.KL = 4;
    radSP.SC = 6;    
    
    multiplier = 1.0;
    
elseif strcmp(op_amp,'1.5x')
    %load image of KLactis cell - this was made in BMH_20130916_analysis_klac.m
    load([ipdir, 'circKL_15x.mat'])
    circSP.KL = circ;
    %load image of S.Cerevisiae cell (made by Jacob)
    load([ipdir,'circSC_15x.mat'])
    circSP.SC = circ;
    
    %Size of image of individual cell to analyze
    sizSP.KL = [18,18];
    sizSP.SC = [25,25];
    %radius of cell for nuclear enrichment calculation
    radSP.KL = 8;
    radSP.SC = 11;
    multiplier = 1.5;
else
    'Error - incorect optical amplification parameter'
end

%Note: these are currently hard coded in
%max_shift_rad = 4 * multiplier;
%min_dist_thresh = 3 * multiplier;

%Tracking parameters
%estimated maximum number of tracks
maxdisp = maxdisp_1x * multiplier;

Nchan = length(channels);



siz = sizSP.(species);
rad = radSP.(species);
circ = circSP.(species);
std_thresh = std_threshSP.(species);
positions = length(pos_fnames);

if length(channel_to_image) == 2
    channels_to_image = channel_to_image;

    for nn = 1:positions
        pos_fnames_nn = pos_fnames{nn};
        for ch = 1:length(channels_to_image)
            channel_to_image = channels_to_image{ch};
            %get filenames and put them into a strucure
            images0 = get_image_fnames(fname_conv,imdir,channel_to_image,pos_fnames_nn,Nchan);
            
            images.(channel_to_image) = images0.(channel_to_image);
            %get times
            [time_inds, time_valsCH.(channel_to_image)] = get_image_times(fname_conv,imdir,channel_to_image,pos_fnames_nn,images,time_calc);
            %time_inds should be the same - time vals might be slightly
            %different
        end

        %Make average time vals? 
        time_vals = zeros(size(time_valsCH.(channels_to_image{1})));
        for ch = 1:length(channels_to_image)
            time_vals = time_vals + time_valsCH.(channels_to_image{ch});
        end
        time_vals = time_vals/length(channels_to_image);
        
        %get tracks
        tracks = cell_track_data(imdir,images,time_inds,time_vals,channels_to_image,circ,imbg,siz,storeim,rad,std_thresh,maxdisp);

        %add position information for each track
        [tracks.pos] = deal(nn);

        %store data
        tracks_vec(nn).tracks = tracks;

        %Starting point if you just want graphs
        %tracks = tracks_vec(nn).tracks;

        plot_each_pos = 0;

        if plot_each_pos == 1
            nTracks = length(tracks);
            figure(1)  
            im = imread([imdir,images(1).(channel_to_image)]); 
            %subplot(2,2,nn)
            hold all
            imagesc(im)
            %imshow(im);

            for jj = 1:nTracks
                x_vec = [tracks(jj).Cxloc];
                y_vec = [tracks(jj).Cyloc];
                %imagesc seems to flip the axes
                plot(y_vec,x_vec,'LineWidth',3)
            end

            [mean_nf, std_nf] = nf_calcs(time_vals,tracks)

            figure(2)
            %subplot(2,2,nn)
            hold all
            for jj = 1:nTracks
                nf_vec = [tracks(jj).nf];
                t_vec = [tracks(jj).times];
                plot(t_vec,nf_vec)
            end




            plot(time_vals,mean_nf, 'k','LineWidth', 3)
            xlabel('time (m)')
            ylabel('Nuclear Localization')
            axis([0,110,1.0,5])


            figure(3)
            %subplot(2,2,nn)
            hold on
            errorbar(time_vals,mean_nf,std_nf)
            xlabel('time (m)')
            ylabel('Nuclear Localization')
            axis([0,50,1.0,5])

            figure(4)
            hold on;
            alpha = 0.2;
            acolor = 'c';
            fill([time_vals fliplr(time_vals)],[mean_nf'+std_nf' fliplr(mean_nf'-std_nf')],acolor, 'FaceAlpha', alpha,'linestyle','none');
            plot(time_vals,mean_nf,acolor,'linewidth',1.5); % change color or linewidth to adjust mean line

        end

    end

end

all_tracks = [];
for nn = 1:length(tracks_vec)
    all_tracks = [all_tracks,tracks_vec(nn).tracks];
end
%right now not averaging time when multiple photos taken - just taking
%last time.  
all_times = time_vals;


%figure(1)
%suptitle('Cell Tracks')     
%figure(2)
%suptitle('Nuclear Localization in response to Glucose loss')        
%figure(3)
%suptitle('Mean Nuclear Localization in response to Glucose loss')
     
%profile off
end

