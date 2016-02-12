function [all_tracks, all_times] = time_series_analysis(analysis_params, usedLED)
% Code for quantifying KL or SC nuclear localization using 1x or 1.5x
% optical zoom.  
% Change jacobs image code output to sort in image order. 
% Too many tracks - something wierd going on
% Had to use ".tiff" for several images
% Note: to make training image use ginput function

%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number, then that is interpreted as a dt between images.  If
%the input is a filename, then that is interpreted as metadata that gives
%exact time values.  



storeim = analysis_params.storeim;
fname_conv = analysis_params.fname_conv;
circ = analysis_params.circ;
siz = analysis_params.siz;
rad = analysis_params.rad;
maxdisp = analysis_params.maxdisp;
track_memory = analysis_params.track_memory;
min_points_for_traj = analysis_params.min_points_for_traj;
%max_shift_rad = analysis_params.max_shift_rad;
%min_dist_thresh = analysis_params.min_dist_thresh;
std_thresh = analysis_params.std_thresh;
thresh = analysis_params.thresh;
pos_fnames = analysis_params.pos_fnames;
channels = analysis_params.channels;
channels_to_image = analysis_params.channels_to_image;
time_calc = analysis_params.time_calc;
imbg = analysis_params.imbg;
maxcells = analysis_params.maxcells;
ipdir = analysis_params.ipdir;
imdir = analysis_params.imdir;

positions = length(pos_fnames);
Nchan = length(channels);

for nn = 1:positions
    thePos = nn % TEST
    pos_fnames_nn = pos_fnames{nn}; % looping through each field of view
    
    if length(channels_to_image) == 2 % if there are 2 channels
        channels_to_image = channels_to_image;
        for ch = 1:length(channels_to_image)
            channel_to_image = channels_to_image{ch};
            %get filenames and put them into a strucure
            images0 = get_image_fnames(fname_conv,imdir,channel_to_image,pos_fnames_nn,Nchan);
            
            images.(channel_to_image) = images0.(channel_to_image);
            %get times
            [time_inds, time_valsCH.(channel_to_image)] = get_image_times(fname_conv,imdir,channel_to_image,pos_fnames_nn,images,time_calc, usedLED);
            %time_inds should be the same - time vals might be slightly
            %different
        end

        %Make average time value for all channels at one position and time
        time_vals = zeros(size(time_valsCH.(channels_to_image{1})));
        for ch = 1:length(channels_to_image)
            time_vals = time_vals + time_valsCH.(channels_to_image{ch});
        end
        time_vals = time_vals/length(channels_to_image);
        channel_to_image = channels_to_image; %convert back to cell so that subsequent subroutines recognize there are two channels. 
        
    elseif length(channels_to_image) == 1 % if there is only 1 channel 
        channel_to_image = channels_to_image{1};       
        %get filenames for 1 FOV and put them into a strucure
        images = get_image_fnames(fname_conv,imdir,channel_to_image,pos_fnames_nn,Nchan); 
    
        %get times
        [time_inds, time_vals] = get_image_times(fname_conv,imdir,channel_to_image,pos_fnames_nn,images,time_calc, usedLED); % this is just for 1 position
    
    end
    time_vals = time_vals/length(channels_to_image);
     

%     else 
%                
%         %get filenames and put them into a strucure
%         images = get_image_fnames(fname_conv,imdir,channel_to_image,pos_fnames_nn,Nchan);
%     
%         %get times
%         [time_inds, time_vals] = get_image_times(fname_conv,imdir,channel_to_image,pos_fnames_nn,images,time_calc);
%     
%     end
    %get tracks

    cell_find_params.circ = circ;
    cell_find_params.imbg = imbg;
    cell_find_params.siz = siz;
    cell_find_params.storeim = storeim;
    cell_find_params.rad = rad;
    cell_find_params.std_thresh = std_thresh;
    cell_find_params.thresh = thresh;
    cell_find_params.coarse_smooth = analysis_params.coarse_smooth;
    cell_find_params.local_smooth = analysis_params.local_smooth; 
    cell_find_params.maxcells = maxcells;
    cell_find_params.deconvlucy_iterations = analysis_params.deconvlucy_iterations;
    cell_find_params.close_max = analysis_params.close_max;
    cell_find_params.far_max = analysis_params.far_max;
    cell_find_params.ne_pixels = analysis_params.ne_pixels;
    cell_find_params.edge_margin = analysis_params.edge_margin;


    timecoursedata = CovMapCells(imdir, images, time_vals, channels_to_image, cell_find_params);
    
    display('finished with CovMapCells');
    
    track_params.time_inds = time_inds;
    track_params.maxdisp = maxdisp;
    track_params.track_memory = track_memory;
    if length(time_inds)<min_points_for_traj
        track_params.min_points_for_traj = 1;
    else
        track_params.min_points_for_traj = min_points_for_traj;
    end
    track_params.maxcells = maxcells;
    track_params.ipdir = ipdir;
    tracks = trackIDL(timecoursedata,track_params);
    %note: tracks use index of times rather than actual times
    %filter const bright cells
    tracks = filterTracks(tracks);
    
    %add position information for each track
    [tracks.pos] = deal(nn);
    
    %store data
    tracks_vec(nn).tracks = tracks;

    display('found tracks')
    
    %Starting point if you just want graphs
    %tracks = tracks_vec(nn).tracks;

%     plot_each_pos = 0;
% 
%     if plot_each_pos == 1
%         nTracks = length(tracks);
%         figure(1)  
%         im = imread([imdir,images(1).(channel_to_image)]); 
%         %subplot(2,2,nn)
%         hold all
%         imagesc(im)
%         %imshow(im);
% 
%         for jj = 1:nTracks
%             x_vec = [tracks(jj).Cxloc];
%             y_vec = [tracks(jj).Cyloc];
%             %imagesc seems to flip the axes
%             plot(y_vec,x_vec,'LineWidth',3)
%         end
% 
%         [mean_nf, std_nf] = nf_calcs(time_vals,tracks)
% 
%         figure(2)
%         %subplot(2,2,nn)
%         hold all
%         for jj = 1:nTracks
%             nf_vec = [tracks(jj).nf];
%             t_vec = [tracks(jj).times];
%             plot(t_vec,nf_vec)
%         end
% 
% 
% 
% 
%         plot(time_vals,mean_nf, 'k','LineWidth', 3)
%         xlabel('time (m)')
%         ylabel('Nuclear Localization')
%         axis([0,110,1.0,5])
% 
% 
%         figure(3)
%         %subplot(2,2,nn)
%         hold on
%         errorbar(time_vals,mean_nf,std_nf)
%         xlabel('time (m)')
%         ylabel('Nuclear Localization')
%         axis([0,50,1.0,5])
% 
%         figure(4)
%         hold on;
%         alpha = 0.2;
%         acolor = 'c';
%         fill([time_vals fliplr(time_vals)],[mean_nf'+std_nf' fliplr(mean_nf'-std_nf')],acolor, 'FaceAlpha', alpha,'linestyle','none');
%         plot(time_vals,mean_nf,acolor,'linewidth',1.5); % change color or linewidth to adjust mean line
% 
%     end

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

