function all_tracks = LH20150527_CRZ1_noCa2_w1488_TIRF_low_s1_t1()
%This is a template for processing nuclear localization data from a set of
%individual .tiff files taken with micromanager.  

%% Inputs

%image processing directory
%ipdir = '/Users/liamholt/Documents/MATLAB/HESlab_nuclear_localization/'
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path 
path(ipdir,path)
analysis_params.ipdir = ipdir;

% Filename convention
% 'JSO','Micromanager', or 'HCS_Nikon', 'Metamorph'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files
fname_conv = 'Metamorph';
analysis_params.fname_conv = fname_conv;

% ?? doublecheck this
analysis_params.storeim = 1;

%Base directory for image files
%base_dir = '/Users/liamholt/Documents/MATLAB/0-CRZ1-test/20150527-CRZ1-noCa2_w1488 TIRF low_'
base_dir = 'C:\Users\Ben\Documents\Data\0-CRZ1-test\20150527-CRZ1-noCa2_w1488 TIRF low_'
%filenames for each phase of imaging.  If only one phase of imaging code
%will still work - this basically is the directory that Micromanager uses
%to store all image files from a certain set of images. 
%example: 
%
%phase = {'Exp'}
%imdirPhase.Exp = [base_dir,'Exp',filesep]
%

phases =  {'Exp'}%,'Post_p2'} %,'Post'} 
%These are the absolute times at which each phase starts.
shift_timing = [0]    

%imdirPhase.Exp = [base_dir,'Exp',filesep]
imdirPhase.Exp = base_dir;

load([ipdir, 'circSC_60x_4x4.mat'])
circ_out = circ;
siz = [16,16];
rad = 8;
std_thresh = 0.2;
maxdisp = 2;

% image of typical cell
analysis_params.circ = circ;

%Size of image of identified cell to extract for further analysis
analysis_params.siz = siz;

%Typical radius of cell. 
analysis_params.rad = rad;

%number of pixels that smoothing uses
analysis_params.coarse_smooth = 25;
analysis_params.local_smooth = 3; 
%number greater than the maximum number of cells in an image.  Used for
%both analysis and tracking
analysis_params.maxcells = 200;

%Number of iterations to use in the deconvlucy algorithm.
analysis_params.deconvlucy_iterations = 5;

%Parameters for finding maxima
analysis_params.close_max = 1;
analysis_params.far_max = 6;

%Number of pixels to use in calculating nuclear enrichment - 
analysis_params.ne_pixels = 5;

%Margin around edge of an image for which cells are removed;
analysis_params.edge_margin = 10;


%Tracking Parameter: estimated maximum displacement
analysis_params.maxdisp = maxdisp;

%Tracking Parameter: 
analysis_params.track_memory = 2;

%Tracking Parameter
analysis_params.min_points_for_traj = 4;

%std thresh for calling a peak to find a cell
analysis_params.std_thresh = std_thresh;

%total number of channels of images taken
analysis_params.channels = {'YFP'};

%Channels we want to process images from.
%If two channels, use a cell.  If just one channel just list it as a text variable i.e. 'RFP'. 
%channel_to_image = 'RFP'
channel_to_image = 'GFP';  
analysis_params.channel_to_image = channel_to_image;

fname_save = 'example_processed_data.mat';


%% timestep method
%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number then that is interpreted as a dt in minutes between images
%If the input is a filename, then that is interpreted as metadata that gives
%exact time values.  

%time_calc_phase.Pre = 5
%time_calc_phase.Post = 4
%time_calc_phase.Pre = [imdirPhase.Pre, 'acqdat.txt']
%time_calc_phase.Post = [imdirPhase.Post, 'acqdat.txt']
%extract actual times from metadata in micromanager images
%addpath('C:/Users/Ben/Documents/GitHub/image_analysis/jsonlab');
%time_calc_phase.Pre  = 'metadata.txt'
%time_calc_phase.Post = 'metadata.txt'
generate_metadata_parsed = 0;
metadata_conv_fname = [ipdir,'times_from_umanager_metadata.py'];
if generate_metadata_parsed ==1;
    for ph = [1:length(phases)]
          phase = phases{ph};
          imdir = imdirPhase.(phase);
          phase = ['python ', metadata_conv_fname, ' ', imdir]
          system(['python ', metadata_conv_fname, ' ', imdir])
    end
end
time_calc_phase.Exp =  .5;
%time_calc_phase.Post_p2 = 'metadata_parsed.txt';





%% Background image

bgimg = 0; 
%1 if you have a background image, 0 if you don't. 
%Gets background image depending on channel to image
%Collect background images using micromanager 
if bgimg == 0  %set default to 1
    if iscell(channel_to_image)    
        ch2i = channel_to_image;
    else
        ch2i = {channel_to_image};
    end
    
    for jj = 1:length(ch2i)
        %ch2i_txt = ch2i{jj}
        imbg.(ch2i{jj}) = 1;
    end
    
else
    %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
    if iscell(channel_to_image)    
        ch2i = channel_to_image;
    else
        ch2i = {channel_to_image};
    end
   
    for jj = 1:length(ch2i)
        ch2i_txt = ch2i{jj}
        if strcmp(ch2i_txt,'RFP')
             imbg_jj = imread([base_dir,'BG',filesep,'img_000000000_RFP_001.tif']);
        elseif strcmp(ch2i_txt,'YFP')
             imbg_jj = imread([base_dir,'BG',filesep,'img_000000000_YFP_001.tif']);
            %imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\YFP_BG\img_000000000_Default_000.tiff');
        else
             'Error: imbg not assigned'
        end
    
        if strcmp(fname_conv,'JSO')
            imbg_jj = imbg_jj';  %micromanager images are the inverse of images collected by JSO image collection scripts.
        end
    
        %Convert imbg to double and median filter
        imbg_jj = double(imbg_jj);
        coarse_smooth = 25;
        %smooths background image using course_smooth parameter.  Boundary
        %conditions are symmetric because default 0 bc's causes strange artifacts
        %on the edges.  For these background images symmetric BCs is a good
        %assumption
        imbg_jj = medfilt2(imbg_jj,[coarse_smooth,coarse_smooth],'symmetric');
        imbg.(ch2i_txt) = imbg_jj;
    end
end

analysis_params.imbg = imbg;


%% Plate information
% Names come from using HCS site generator and micromanager
% all locations had 4 sites
% Glucose Dropout (0.1 % from 2%) followed by Osmo shock (0.35M) at various times.
% osmo balanced with sorbitol relative to osmotic stress condition.
%
%
%A8	t9: GD + sorb   
%B8	t9: GD      
%C8	t9: SDC + sorb
%D8	t9: GD  	t15: sorb
%E8	t9: GD      t21: sorb
%F8	t9: GD      t33: sorb
%G8	t9: GD      t45: sorb
%H8	t9: GD      t69: sorb



%should also be able to delete individual images from a series to remove
%those and code will still run.  


% Using JSO's matlab image code, the folder names are different - here is
% an example of the position vector in that case
 posvec = {'s1','s2';
 's3','s4'};
% 'p11','p12','p13';
% 'p16','p17','p18';
% 'p21','p22','p23';
% 'p26','p27','p28'};

%Remove various sites from list 

%note that the difference between python site numbering in filenames 
% and matlab numbering in the position vector means you must add one to the 
% site number to reference the vector 
%e.g. to remove Remove A1-Site_1
%
%posvecSP.SC{1,2} = 'NA'
%


%Obtain and store data for each dose 

all_tracks_vec = [];
all_times_vec = [];
Nconditions = size(posvec,2);
for jj = 1:Nconditions
   for ph = 1:length(phases);
        phase = phases{ph};
        analysis_params.imdir = imdirPhase.(phase);
        analysis_params.time_calc = time_calc_phase.(phase);
        pos_fnames = posvec(jj,:); 
        %remove bad positions (NAs)
        analysis_params.pos_fnames = pos_fnames(~strcmp(pos_fnames,'NA'));
        %main function for data processing
        [tracks,times] = time_series_analysis(analysis_params);
        all_tracks.(phase) = tracks;
        all_times.(phase) = times + shift_timing(ph);
   end
   all_tracks_vec{jj} = all_tracks;
   all_times_vec{jj} = all_times;
end
%store data
save([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

return
