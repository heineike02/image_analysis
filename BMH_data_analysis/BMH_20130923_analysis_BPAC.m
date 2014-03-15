% Incorporate Cell Tracking
% Try with K. Lactis
% Visualize individual and mean nuclear localization for glucose shift
% before and after Glucose deprivation in S.C. and K. Lactis
% Repeat for BPAC (in another file)

% Had to use ".tiff" for several images
% Note: to make training image use ginput function

profile on

%load image - first image is S. Cerevisiae pre-induction
%Windows: 
ipdir = 'C:\Users\Ben\Google Drive\HES_Lab_Share\Scripts\JSO_Image_Analysis\BMH_Image_Analysis\'
imdirSP.KL = 'C:\Users\Ben\Documents\Data\PKA_Project\20130923\Blue_light_Strain_KL_BPAC\'
imdirSP.SC = 'C:\Users\Ben\Documents\Data\PKA_Project\20130923\Blue_light_Strain_SC_BPAC\'

%Gets background image for RFP
imbg = imread('C:\Users\Ben\Documents\Data\PKA_Project\20130923\RFP_BG100ms\img_000000000_Default_000.tif');
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



%Mac
%ex: imdir = '/Users/ucsf/Google Drive/UCSF/ElSamad_Lab/PKA/WetLab/20130716/Gluc_pre_ind/Pos1/'
storeim = 1;
%load image of KLactis cell - this was made in BMH_20130916_analysis_klac.m
load([ipdir, 'circKL.mat'])
circSP.KL = circKL;

%load image of S.Cerevisiae cell (made by Jacob)
load([ipdir(1:end-19),'circ.mat'])
circSP.SC = circ;

%[celldata,N] = FindCellsBMH(im, circ, bf_mask, im_bg , siz, storeim)

%Size of image of individual cell to analyze
sizSP.KL = [15,15];
sizSP.SC = [17,17];
%radius of cell for nuclear enrichment calculation
radSP.KL = 4;
radSP.SC = 6;
%std thresh for calling a peak to find a cell
%0.16 seemed to work best for k.lactis on 16JUL images
std_threshSP.KL = 0.16;
std_threshSP.SC = 0.2;

max_shift_rad = 4;
min_dist_thresh = 3;

%estimated maximum number of tracks
max_tracks = 100;


%There are four positions per condition.  

species_cell = {'SC'}
channels = {'RFP', 'YFP'}
%Eventually should track cells using both channels, but for now RFP is the
%best channel
cell_track_chan_ind = 1
cell_track_chan = {'RFP'}
Nchan = length(channels)
positions = 4
dt = .5

%initialize data structure to store all data
%tracks_init = struct('Cxloc',zeros(nTimes,1),'Cyloc', zeros(nTimes,1), 'nf', zeros(nTimes,1), 'nmi', zeros(nTimes,1), 'times',zeros(nTimes,1), 'length',0 );
%tracks(1:max_tracks) = tracks_init;
%tracks_vec = struct('tracks',tracks_init);
%tracks_vec = struct('tracks',[]);
%tracks_vec(positions).tracks = []


for sp = 1:length(species_cell)
    species = species_cell{sp};
    siz = sizSP.(species);
    rad = radSP.(species);
    circ = circSP.(species);
    std_thresh = std_threshSP.(species);
    imdir = imdirSP.(species);
    D = dir(imdir);
    %subtract 6 for the non-image files in the directory.
    nTimes = (length(D)-6)/(Nchan*positions);
    times = dt*[0:nTimes-1];


        
    %initializing plotting variables
    nf_vec_mean = zeros([nTimes,positions]);
    
    for nn = 1:positions
        
        pos = num2str(nn);
                
        for mm = 1:Nchan
           images.(channels{mm}) = cell([nTimes,1]);        
           for jj = 1:length(D)
               fname = D(jj).name;
               %check correct channel and position
               if strfind(fname,[channels{mm},'_p',pos]) > 0
                  %find time point and convert to index.  Channel ordering has to be in correct order (channel one t1, channel two t2 etc: 
                  t_ind = regexp(fname, 't(?<time>\d+)','names');
                  t_ind = str2double(t_ind.time);
                  kk = (t_ind + Nchan-mm)/Nchan;
                  images.(channels{mm}){kk} = D(jj).name;
               end
           end
        end
        
        
        %Get data for each position
        timecoursedata = CovMapCellsBMH(imdir ,images, channels, circ, imbg , siz , storeim, rad, std_thresh, times);
        %Just using RFP data to track cells
        timecoursedata_RFP = timecoursedata(:,cell_track_chan_ind);
       
        maxdisp = 4;
        tracks = trackIDL_BMH(timecoursedata_RFP,times, maxdisp);
              
        %filter const bright cells
        tracks = filterTracksBMH(tracks);
        %store data
        tracks_vec(nn).tracks = tracks;
        
        %Starting point if you just want graphs
        %tracks = tracks_vec(nn).tracks;
        
        nTracks = length(tracks);
        figure(1)
        im = imread([imdir,images.(cell_track_chan{1}){1}]); 
        subplot(2,2,nn)
        hold all
        imagesc(im)
        
        
        for jj = 1:nTracks
            x_vec = [tracks(jj).Cxloc];
            y_vec = [tracks(jj).Cyloc];
            %imagesc seems to flip the axes
            plot(y_vec,x_vec)
        end
        

        
        figure(2)
        subplot(2,2,nn)
        hold all
        for jj = 1:nTracks
            nf_vec = [tracks(jj).nf];
            t_vec = [tracks(jj).times];
            plot(t_vec,nf_vec)
        end
        
        
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
        
        plot(times,mean_nf, 'k','LineWidth', 3)
        xlabel('time (m)')
        ylabel('Nuclear Localization')
        axis([0,40,1.0,4])
        
        
        figure(3)
        subplot(2,2,nn)
        errorbar(times,mean_nf,std_nf)
        xlabel('time (m)')
        ylabel('Nuclear Localization')
        axis([0,40,1.0,4])
        
               

        
    end

      
end

figure(1)
suptitle('Cell Tracks')     
figure(2)
suptitle('Nuclear Localization: Light on between 4-7min')        
figure(3)
suptitle('Mean Nuclear Localization: Light on between 4-7min')
        

profile off