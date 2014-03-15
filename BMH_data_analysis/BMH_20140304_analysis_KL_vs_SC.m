function BMH_20140304_analysis_KL_vs_SC()
% Glucose Dropout and Oxidative stress in Lactis and Cerevisiae
% Visualize individual and mean nuclear localization each condition
% in S.C. and K. Lactis

% Used 1.0x optical zoom.  
% Need to have argument for channel to measure / output
% Change jacobs image code output to sort in image order. 
% Too many tracks - something wierd going on

% Had to use ".tiff" for several images
% Note: to make training image use ginput function

profile on

%load image - first image is S. Cerevisiae pre-induction
%Windows: 
ipdir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\'

%Mac
%ex: imdir = '/Users/ucsf/Google Drive/UCSF/ElSamad_Lab/PKA/WetLab/20130716/Gluc_pre_ind/Pos1/'
storeim = 1;

% Filename convention
% 'JSO' or 'Micromanager'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files
fname_conv = 'JSO'

% optical amplification 
% 1x or 1.5x
op_amp =  '1.5x'
storeim = 1;
imbg = 1;

%[celldata,N] = FindCellsBMH(im, circ, bf_mask, im_bg , siz, storeim)

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

%[celldata,N] = FindCellsBMH(im, circ, bf_mask, im_bg , siz, storeim)

%std thresh for calling a peak to find a cell
%0.16 seemed to work best for k.lactis on 16JUL images
std_threshSP.KL = 0.16;
std_threshSP.SC = 0.2;

max_shift_rad = 4 * multiplier;
min_dist_thresh = 3 * multiplier;

%Tracking parameters
%estimated maximum number of tracks
max_tracks = 100;
maxdisp = 4 * multiplier;

%input information from experiment here
species_cell = {'KL'} %{'SC'}  %
phases = {'Pre','Post'}
shift_spacing = [0,4]

%imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_pre_1\'
%imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\'

base_dir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140304\'
imdirPhase.Pre = [base_dir,'Ox_Gluc_Drop_Pre\']
imdirPhase.Post = [base_dir,'Ox_Gluc_Drop_Post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc
%pos_fnamesSP.SC = {'Pos0','Pos1','Pos2','Pos3','Pos4','Pos5','Pos6'};
%pos_fnamesSP.KL = {'Pos3','Pos4','Pos5'};
pos_fnamesSP.SC = {'p1','p2','p3'}
pos_fnamesSP.KL = {'p10','p11','p12'} 
channels = {'YFP','RFP','BF'}
channel_to_image = 'YFP'
Nchan = length(channels)
dt = 1.0

%Gets background image depending on channel to image
%Collect background images using micromanager 
%imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
if strcmp(channel_to_image,'RFP')
    imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140304\Ox_Gluc_Drop_Pre\RFP_p12_t11.tiff');
elseif strcmp(channel_to_image,'YFP')
    imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140304\Ox_Gluc_Drop_Pre\YFP_p1_t3.tiff');
else
    'Error: imbg not assigned'
end
%Default imbg is 1
%imbg = 1;
%Convert imbg to double and median filter
imbg = double(imbg);
coarse_smooth = 25;
%smooths background image using course_smooth parameter.  Boundary
%conditions are symmetric because default 0 bc's causes strange artifacts
%on the edges.  For these background images symmetric BCs is a good
%assumption
imbg = medfilt2(imbg,[coarse_smooth,coarse_smooth],'symmetric');



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
                    acolor = 'b';
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
            
            figure(6)
            hold on;
            alpha = 0.2;
            acolor = 'b';
            fill([times fliplr(times)],[mean_nf'+std_nf' fliplr(mean_nf'-std_nf')],acolor, 'FaceAlpha', alpha,'linestyle','none');
            plot(times,mean_nf,acolor,'linewidth',1.5); % change color or linewidth to adjust mean line
            xlabel('time(min)')
            ylabel('Mean Nuclear Localization')
            title('K.Lac Ox Stress')
            
            figure(7)
            hold all
            for nn = 1:length(all_tracks)
               nf_vec = [all_tracks(nn).nf];
               t_vec = [all_tracks(nn).times];
               plot(t_vec,nf_vec)
            end
            xlabel('time(min)')
            ylabel('Nuclear Localization')
            title('K.Lac Ox Stress')
            
        end


    end


    figure(1)
    suptitle('Cell Tracks')     
    figure(2)
    suptitle('Nuclear Localization in response to Glucose loss/addition')        
    figure(3)
    suptitle('Mean Nuclear Localization in response to Glucose loss/addition')
        

elseif strcmp(fname_conv,'Micromanager')


%See earlier date files to see how to manage data with Micromanager file
%naming conventions

profile off


end

end


function [tracks,times,images] = cell_track_data(imdir,channel_to_image,pos_fnames_nn,Nchan,shift_spacing_ph,dt,circ,imbg,siz,storeim,rad,std_thresh,maxdisp,t);

    images = struct();
    files=dir(char(strcat(imdir,channel_to_image,'_',pos_fnames_nn,'_t*.tiff'))); 
    nTimes = size(files,1);
    %Get initial number
    init_fname = files(1).name;
    init_ind = strfind(init_fname,'_t');
    dot_ind = strfind(init_fname,'.');
    init_num = str2num(init_fname(init_ind+2:(dot_ind-1)));
    %picks out the first number and checks what it is modulo the number of channels 
    channel_mod = mod(init_num,Nchan)+Nchan*(mod(init_num,Nchan)==0) ;
    %This is a pain because the numbering is not 01,02,03 and
    %messes up the ordering of the files. 
    for mm = 1:nTimes 
       images(mm).(channel_to_image) = char(strcat(channel_to_image,'_',pos_fnames_nn,'_t',int2str(channel_mod+(mm-1)*Nchan),'.tiff'));
    end
    
    %Shift times for each phase
    t = t + shift_spacing_ph;
    times = dt*[0:nTimes-1]+t;
    t = times(end);
    
    timecoursedata = CovMapCellsBMH(imdir ,images, {channel_to_image}, circ, imbg , siz , storeim, rad, std_thresh, times);
    tracks = trackIDL_BMH(timecoursedata,times, maxdisp);

    %filter const bright cells
    tracks = filterTracksBMH(tracks);

end


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



