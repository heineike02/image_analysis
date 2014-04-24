function all_tracks = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species_cell,phases, shift_spacing,imdirPhase, maxdisp_1x,pos_fnamesSP,channels,channel_to_image,dt,imbg,summary_title,acolor_summary)
% Code for quantifying KL or SC nuclear localization using 1x or 1.5x
% optical zoom.  
% Change jacobs image code output to sort in image order. 
% Too many tracks - something wierd going on
% Had to use ".tiff" for several images
% Note: to make training image use ginput function

profile on
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
    load([ipdir(1:end-19),'circSC_1x.mat'])
    circSP.SC = circ;
    
    %Size of image of individual cell to analyze
    sizSP.KL = [15,15];
    sizSP.SC = [17,17];
    %radius of cell for nuclear enrichment calculation
    radSP.KL = 4;
    radSP.SC = 6;    
    
    multiplier = 1.0
    
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

Nchan = length(channels)

if strcmp(fname_conv,'JSO')

    
    for sp = 1:length(species_cell)
        species = species_cell{sp};
        siz = sizSP.(species);
        rad = radSP.(species);
        circ = circSP.(species);
        std_thresh = std_threshSP.(species);
        pos_fnames = pos_fnamesSP.(species);
        positions = length(pos_fnames)
        
        t = 0
        for ph = 1:length(phases);
            phase = phases{ph};
            shift_spacing_ph = shift_spacing(ph);
            imdir = imdirPhase.(phase);
            
            for nn = 1:positions
                pos_fnames_nn = pos_fnames(nn);
                [tracks, times, images] = cell_track_data(imdir,channel_to_image,pos_fnames_nn,Nchan,shift_spacing_ph,dt,circ,imbg,siz,storeim,rad,std_thresh,maxdisp,t);
                nTimes = length(times);
                
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

                    [mean_nf, std_nf] = nf_calcs(times,tracks)

                    figure(2)
                    %subplot(2,2,nn)
                    hold all
                    for jj = 1:nTracks
                        nf_vec = [tracks(jj).nf];
                        t_vec = [tracks(jj).times];
                        plot(t_vec,nf_vec)
                    end

                    
                    

                    plot(times,mean_nf, 'k','LineWidth', 3)
                    xlabel('time (m)')
                    ylabel('Nuclear Localization')
                    axis([0,110,1.0,5])


                    figure(3)
                    %subplot(2,2,nn)
                    hold on
                    errorbar(times,mean_nf,std_nf)
                    xlabel('time (m)')
                    ylabel('Nuclear Localization')
                    axis([0,50,1.0,5])

                    figure(4)
                    hold on;
                    alpha = 0.2;
                    acolor = 'c';
                    fill([times fliplr(times)],[mean_nf'+std_nf' fliplr(mean_nf'-std_nf')],acolor, 'FaceAlpha', alpha,'linestyle','none');
                    plot(times,mean_nf,acolor,'linewidth',1.5); % change color or linewidth to adjust mean line

                end
                
                
                
            end
            
            t = times(end);
            
            all_tracks = [];
            for nn = 1:length(tracks_vec)
                all_tracks = [all_tracks,tracks_vec(nn).tracks];
            end
            
            [mean_nf, std_nf] = nf_calcs(times,all_tracks);
            
            figure(9)
            hold on;
            alpha = 0.2;
            fill([times fliplr(times)],[mean_nf'+std_nf' fliplr(mean_nf'-std_nf')],acolor_summary, 'FaceAlpha', alpha,'linestyle','none');
            plot(times,mean_nf,acolor_summary,'linewidth',1.5); % change color or linewidth to adjust mean line
            xlabel('time(min)')
            ylabel('Mean Nuclear Localization')
            title(summary_title)
            
            figure(7)
            hold all
            for nn = 1:length(all_tracks)
               nf_vec = [all_tracks(nn).nf];
               t_vec = [all_tracks(nn).times];
               plot(t_vec,nf_vec)
            end
            xlabel('time(min)')
            ylabel('Nuclear Localization')
            title(summary_title)
            
        end


    end


    figure(1)
    suptitle('Cell Tracks')     
    figure(2)
    suptitle('Nuclear Localization in response to Glucose loss')        
    figure(3)
    suptitle('Mean Nuclear Localization in response to Glucose loss')
        

elseif strcmp(fname_conv,'Micromanager')


%See earlier date files to see how to manage data with Micromanager file
%naming conventions

profile off


end

end

