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
imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20130923\Gluc_Drop_Post1_Pre2_1\'
imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20130923\Gluc_Drop_Post2_1\'

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



%Note: at Image21,everything starts to go out of focus.  Should only use up
%to 20


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

species_cell = {'SC'}  %,'KL'}
phases = {'Pre','Post'}
pos_fnamesSP.SC = {'Pos0'} %, 'Pos1','Pos2'};
%pos_fnamesSP.KL = {'Pos3','Pos4','Pos5'};
pos_fnamesSP.KL = {'Pos4'};
channels = {'RFP_001', 'BF_001','BF_000'}
%Eventually should track cells using both channels, but for now RFP is the
%best channel
cell_track_chan_ind = 1
Nchan = length(channels)
positions = 1 %3
dt = 1.0

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
    pos_fnames = pos_fnamesSP.(species);
   
    for ph = 1:length(phases);
        phase = phases{ph};
        
        for nn = 1:positions
            imdir = [imdirPhase.(phase),pos_fnames{nn},'\'];
            D = dir(imdir);
            

            
            %remove last 9 points from the Post images because they went out of focus
            %Shift times for Post addition of glucose to start at 7:30
            if strcmp(phase,'Post')==1
                throwout = 9;
                shift = 7.5;
            else 
                throwout = 0;
                shift = 0;
            end
            
            %subtract 3 for the non-image files listed in D.
            nTimes = (length(D)-3-throwout*Nchan)/(Nchan);
            times = dt*[0:nTimes-1]+shift;

            %initializing plotting variables
            nf_vec_mean = zeros([nTimes,positions]);
            
            for mm = 1:Nchan
               images.(channels{mm}) = cell([nTimes,1]);        
               kk = 1;
               for jj = 1:length(D);
                   fname = D(jj).name;
                   %Note: only works because images are in order of time
                   if strfind(fname,[channels{mm}]) > 0
                      images.(channels{mm}){kk} = D(jj).name;
                      kk = kk + 1;
                   end
               end
            end

            %Get data for each position
            %For now only using RFP images, not using BF or out of focus BF
            %throwout images for Post files due to drifting out of focus: 
            images.RFP = {images.RFP_001{1:(end-throwout)}}';

            
            timecoursedata = CovMapCellsBMH(imdir ,images, {'RFP'}, circ, imbg , siz , storeim, rad, std_thresh, times);
            
            maxdisp = 4;
            tracks = trackIDL_BMH(timecoursedata,times, maxdisp);

            %filter const bright cells
            tracks = filterTracksBMH(tracks);
            %store data
            tracks_vec(nn).tracks = tracks;
         

            %Starting point if you just want graphs
            %tracks = tracks_vec(nn).tracks;

            nTracks = length(tracks);
            figure(1)
            im = imread([imdir,images.RFP{1}]); 
            %subplot(3,1,nn)
            hold all
            imagesc(im)
            %imshow(im);

            for jj = 1:nTracks
                x_vec = [tracks(jj).Cxloc];
                y_vec = [tracks(jj).Cyloc];
                %imagesc seems to flip the axes
                plot(y_vec,x_vec,'LineWidth',3)
            end



            figure(2)
            %subplot(3,1,nn)
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
            axis([0,30,1.0,9.0])


            figure(3)
            %subplot(3,1,nn)
            hold on
            errorbar(times,mean_nf,std_nf)
            xlabel('time (m)')
            ylabel('Nuclear Localization')
            axis([0,30,1,7.0])




        end
    end

      
end

figure(1)
suptitle('Cell Tracks')     
figure(2)
suptitle('Nuclear Localization in response to Glucose Dropout')        
figure(3)
suptitle('Mean Nuclear Localization in response to Glucose Dropout')
        

profile off