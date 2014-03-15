% Glucose Dropout in Lactis and Cerevisiae
% 04FEB14
% Visualize individual and mean nuclear localization for glucose shift
% before and after Glucose deprivation in S.C. and K. Lactis

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
species_cell = {'SC'} %{'SC'}  %
phases = {'Pre','Post','Rep'}
shift_spacing = [0,2,3]

%imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_pre_1\'
%imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\'

base_dir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140204\'
imdirPhase.Pre = [base_dir,'Gluc_Drop_Pre\']
imdirPhase.Post = [base_dir,'Gluc_Drop_Post\']
imdirPhase.Rep_ot = [base_dir,'Gluc_Rep_Post_optest\']
imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc
%pos_fnamesSP.SC = {'Pos0','Pos1','Pos2','Pos3','Pos4','Pos5','Pos6'};
%pos_fnamesSP.KL = {'Pos3','Pos4','Pos5'};
pos_fnamesSP.SC = {'p1'}
pos_fnamesSP.KL = {'p4'} 
channels = {'YFP','RFP'}
channel_to_image = 'RFP'
Nchan = length(channels)
dt = 1.0

%Gets background image depending on channel to image
%Collect background images using micromanager 
imbg = imread([base_dir,'YFP_RFP_BG\BG_1\img_000000000_',channel_to_image,'_000.tif']);
%imbg = imread('C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\Pos1\img_000000000_YFP_000.tif');
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

            for nn = 1:positions
                images = struct();
                imdir = imdirPhase.(phase);
                files=dir(char(strcat(imdir,channel_to_image,'_',pos_fnames(nn),'_t*.tiff'))); 
                nTimes = size(files,1);
                %Get initial number
                init_fname = files(1).name;
                init_ind = strfind(init_fname,'_t');
                dot_ind = strfind(init_fname,'.');
                init_num = str2num(init_fname(init_ind+2:(dot_ind-1)));
                %picks out the first number and checks if it is even or odd
                %to determine if the initial time is 1 or 2.
                one_or_two = 2-mod(init_num,2);
                %This is a pain because the numbering is not 01,02,03 and
                %messes up the ordering of the files. 
                for mm = 1:nTimes 
                   images(mm).(channel_to_image) = char(strcat(channel_to_image,'_',pos_fnames(nn),'_t',int2str(one_or_two+(mm-1)*2),'.tiff'));
                end
                
                %Shift times for each phase
                t = t + shift_spacing(ph);
                times = dt*[0:nTimes-1]+t;
                t = times(end);

                timecoursedata = CovMapCellsBMH(imdir ,images, {channel_to_image}, circ, imbg , siz , storeim, rad, std_thresh, times);
                tracks = trackIDL_BMH(timecoursedata,times, maxdisp);

                %filter const bright cells
                tracks = filterTracksBMH(tracks);
                %store data
                tracks_vec(nn).tracks = tracks;


                %Starting point if you just want graphs
                %tracks = tracks_vec(nn).tracks;

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



                figure(2)
                %subplot(2,2,nn)
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
                axis([0,110,1.0,5])


                figure(3)
                %subplot(2,2,nn)
                hold on
                errorbar(times,mean_nf,std_nf)
                xlabel('time (m)')
                ylabel('Nuclear Localization')
                axis([0,50,1.0,5])




            end
        end


    end


    figure(1)
    suptitle('Cell Tracks')     
    figure(2)
    suptitle('Nuclear Localization in response to Glucose loss/addition')        
    figure(3)
    suptitle('Mean Nuclear Localization in response to Glucose loss/addition')
        

elseif strcmp(fname_conv,'Micromanager')



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
    positions = length(pos_fnames)
   
    for ph = 1:length(phases);
        phase = phases{ph};
        
        for nn = 1:positions
            imdir = [imdirPhase.(phase),pos_fnames{nn},'\'];
            D = dir(imdir);
            

            
            %remove last 9 points from the Post images because they went out of focus
            %Shift times for Post addition of glucose to start at 7:30
            if strcmp(phase,'Post')==1
                throwout = 0;
                shift = 4;
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

            % CHANGE:
            %Get data for each position
            %For now only using YFP images, not using BF or out of focus BF
            %throwout images for Post files due to drifting out of focus: 
            images.RFP = {images.RFP{1:(end-throwout)}}';

            % CHANGE:
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
            im = imread([imdir,images.YFP{1}]); 
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



            figure(2)
            %subplot(2,2,nn)
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
            axis([0,50,1.0,3.5])


            figure(3)
            %subplot(2,2,nn)
            hold on
            errorbar(times,mean_nf,std_nf)
            xlabel('time (m)')
            ylabel('Nuclear Localization')
            axis([0,50,1.0,3.5])




        end
    end

      
end


figure(1)
suptitle('Cell Tracks')     
figure(2)
suptitle('Nuclear Localization in response to Glucose loss/addition')        
figure(3)
suptitle('Mean Nuclear Localization in response to Glucose loss/addition')
        


else

    'Error - no filename convention indicated'
    
end


profile off
