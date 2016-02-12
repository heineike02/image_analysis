function all_tracks = SYC_image_analysis()
%This is a template for processing nuclear localization data from a set of
%individual .tiff files taken with micromanager.  

%% Inputs
% Stationary OD
%image processing directory
ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
%adds image analysis directory to path 
path(ipdir,path);
analysis_params.ipdir = ipdir;

% Filename convention
% 'JSO','Micromanager', or 'HCS_Nikon'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files
fname_conv = 'JSO';
analysis_params.fname_conv = fname_conv;

% ?? doublecheck this
analysis_params.storeim = 1;

%%%%%%
%Base directory for image files
base_dir = '\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20150916_light-mediated_shuttling_NLSvar3-29\';
%base_dir = '\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20150828_light-mediated_shuttling_pSYC83_pSYC84\';
%filenames for each phase of imaging.  If only one phase of imaging code
%will still work - this basically is the directory that Micromanager uses
%to store all image files from a certain set of images. 
%example: 
%
%phase = {'Exp'}
%imdirPhase.Exp = [base_dir,'Exp',filesep]
%

%phases = {'Pre'};
phases =  {'Original','Three','Five','Eight','Nine','Eleven','Fourteen', 'Fifteen','Twenty','TwentyFour', 'TwentySeven','TwentyNine'};%{'Post'};
%These are the absolute times at which each phase starts.
% subfolder = {'pSYC83_10minBL_255auInten10_1point5x','pSYC83_nobluelight_1point5x', 'pSYC84_10minBL_255auInten4_OFF_1point5x',...
%     'pSYC84_10minBL_255auInten10_1point5x','pSYC84_nobluelight_1point5x'};
subfolder = {'1','2','3','4','5','6','7','8','9','10','11','12'};
shift_timing = [0 0 0 0 0 0 0 0 0 0 0 0]; %[0,10]; % assuming minutes    

%imdirPhase.Pre = [base_dir,'before_phospH',filesep];
%imdirPhase.Post = [base_dir,'5minBL_1fov_1point5x',filesep] %'\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20150710_PhosphateDeplet_ASOE_TRpanel\phosphateDeplet'
imdirPhase.Original = [base_dir,subfolder{1}, filesep]
imdirPhase.Three = [base_dir, subfolder{2}, filesep]
imdirPhase.Five = [base_dir, subfolder{3}, filesep]
imdirPhase.Eight = [base_dir, subfolder{4}, filesep]
imdirPhase.Nine = [base_dir, subfolder{5}, filesep]
imdirPhase.Eleven = [base_dir, subfolder{6}, filesep]
imdirPhase.Fourteen = [base_dir, subfolder{7}, filesep]
imdirPhase.Fifteen = [base_dir, subfolder{8}, filesep]
imdirPhase.Twenty = [base_dir, subfolder{9}, filesep]
imdirPhase.TwentyFour = [base_dir, subfolder{10}, filesep]
imdirPhase.TwentySeven = [base_dir, subfolder{11}, filesep]
imdirPhase.TwentyNine = [base_dir, subfolder{12}, filesep]
%imdirPhase.Post = [base_dir,'Post',filesep]

% optical amplification 
% 1x or 1.5x
op_amp = '1p5x';
% species 
% S. Cerevisiae 'SC'
% K. Lactis 'KL'
species = 'SC';

%estimated maximum displacement in pixels between frames for cells in a 40x air
%objective at 1x optical amplification. 
maxdisp_1x = 4;

[circ, siz, rad, maxdisp, std_thresh] = species_magnification_params_syc(species, op_amp, ipdir, maxdisp_1x);

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
analysis_params.far_max = 10;

%Number of pixels to use in calculating nuclear enrichment - 
analysis_params.ne_pixels = 5;

%Margin around edge of an image for which cells are removed;
analysis_params.edge_margin = 10;


%Tracking Parameter: estimated maximum displacement
analysis_params.maxdisp = maxdisp;

%Tracking Parameter: 
analysis_params.track_memory = 2;

%Tracking Parameter
analysis_params.min_points_for_traj = 12; %4;

%std thresh for calling a peak to find a cell
analysis_params.std_thresh = std_thresh;

%thresh for identifying cell locations after deconvlucy
analysis_params.thresh = 2; % standard for S. cerevisiae is 2

%total number of channels of images taken
analysis_params.channels = {'RFP'};

%Channels we want to process images from.
%If two channels, use a cell.  If just one channel just list it as a text variable i.e. 'RFP'. 
%channel_to_image = 'RFP'
channels_to_image = {'RFP'};  
analysis_params.channels_to_image = channels_to_image;

fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150914_light-mediated_shuttling_NLSvar3-29\OrigNLSvar3-29_lightmediatedshuttling_20150916_satOD.mat';

%% timestep method
%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number then that is interpreted as a dt in minutes between images
%If the input is a filename, then that is interpreted as metadata that gives
%exact time values.  

% this is the fixed dT option
%time_calc_phase.Pre = 5
%time_calc_phase.Post = 4

%time_calc_phase.Pre = [imdirPhase.Pre, 'acqdat.txt']
%time_calc_phase.Post = [imdirPhase.Post, 'acqdat.txt']

% extract actual times from metadata in Micromanager images
%addpath('C:/Users/Ben/Documents/GitHub/image_analysis/jsonlab');
%time_calc_phase.Pre  = 'metadata.txt'
%time_calc_phase.Post = 'metadata.txt'
% generate_metadata_parsed = 0;
% metadata_conv_fname = [ipdir,'times_from_umanager_metadata.py'];
% if generate_metadata_parsed ==1;
%     for ph = [1:length(phases)]
%           phase = phases{ph};
%           imdir = imdirPhase.(phase);
%           phase = ['python ', metadata_conv_fname, ' ', imdir]
%           system(['python ', metadata_conv_fname, ' ', imdir])
%     end
% end
% time_calc_phase.Pre =  'metadata_parsed.txt';
% time_calc_phase.Post = 'metadata_parsed.txt';
%time_calc_phase.Post_p2 = 'metadata_parsed.txt';

% Is LED used or not - different columns for the metadata
%usedLED = [1 0 1 1 0];
usedLED = [1 1 1 1 1 1 1 1 1 1 1 1];
% this is the JSO version option
%time_calc_phase.Pre = '\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20150716_PhosphateDeplet_pH_ASOE_5TRs\before_phospH\acqdat.txt';
%time_calc_phase.Post = '\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20150830_light-mediated_shuttling_pSYC84\5minBL_1fov_1point5x\acqdat.txt';
time_calc_phase.Original = strcat(base_dir, subfolder{1},'\acqdat.txt');
time_calc_phase.Three = strcat(base_dir, subfolder{2}, '\acqdat.txt');
time_calc_phase.Five = strcat(base_dir, subfolder{3}, '\acqdat.txt');
time_calc_phase.Eight = strcat(base_dir, subfolder{4}, '\acqdat.txt');
time_calc_phase.Nine = strcat(base_dir, subfolder{5}, '\acqdat.txt');
time_calc_phase.Eleven = strcat(base_dir, subfolder{6}, '\acqdat.txt');
time_calc_phase.Fourteen = strcat(base_dir, subfolder{7}, '\acqdat.txt');
time_calc_phase.Fifteen = strcat(base_dir, subfolder{8}, '\acqdat.txt');
time_calc_phase.Twenty = strcat(base_dir, subfolder{9}, '\acqdat.txt');
time_calc_phase.TwentyFour = strcat(base_dir, subfolder{10}, '\acqdat.txt');
time_calc_phase.TwentySeven = strcat(base_dir, subfolder{11}, '\acqdat.txt');
time_calc_phase.TwentyNine = strcat(base_dir, subfolder{12}, '\acqdat.txt');

%% Background image

bgimg = 0; 
%1 if you have a background image, 0 if you don't. 
%Gets background image depending on channel to image
%Collect background images using micromanager 
if bgimg == 0  %set default to 1
    if iscell(channels_to_image)    
        ch2i = channels_to_image;
    else
        ch2i = {channels_to_image};
    end
    
    for jj = 1:length(ch2i)
        %ch2i_txt = ch2i{jj}
        imbg.(ch2i{jj}) = 1;
    end
    
else
    %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
    if iscell(channels_to_image)    
        ch2i = channels_to_image;
    else
        ch2i = {channels_to_image};
    end
   
    for jj = 1:length(ch2i)
        ch2i_txt = ch2i{jj};
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

% This is the Micromanager position indices
% wellvec = {'A9','B9'} %,'C9','D9','E9','F9','G9','H9'};
% Nsites = 4;
% Nwells = length(wellvec);
% for jj = 1:length(wellvec);
%     for kk = 1:Nsites;
%         posvec{jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
%     end
% end

%Remove various sites from list 

%note that the difference between python site numbering in filenames 
% and matlab numbering in the position vector means you must add one to the 
% site number to reference the vector 
%e.g. to remove Remove A1-Site_1
%
%posvecSP.SC{1,2} = 'NA'
%

%should also be able to delete individual images from a series to remove
%those and code will still run.  


% Using JSO's matlab image code, the folder names are different - here is
% an example of the position vector in that case
%posvec = {'p1'};
posvec = {'p1'};
%{'p1','p2','p3','p4';
%     'p5','p6','p7','p8';
%     'p9','p10','p11','p12';
%     'p13','p14','p15','p16';
%     'p17','p18', 'p19', 'p20';
%     'p21','p22','p23', 'p24';
%     'p25','p26','p27', 'p28';
%     'p29','p30','p31','p32'};

%remember that well p19 is bad

Nwells = size(posvec,1);

% to remove just one position - do this, but this doesn't work if you
% remove the whole row
%posvec{5,3} = 'NA';

%Obtain and store data for each dose 

all_tracks_vec = [];
all_times_vec = [];
posvec;
for jj = 1:Nwells
   Wells = jj;
   for ph = 1:length(phases);
        ph
        phase = phases{ph};
        analysis_params.imdir = imdirPhase.(phase);
        analysis_params.time_calc = time_calc_phase.(phase);
        pos_fnames = posvec(jj,:); 
        %remove bad positions (NAs)
        analysis_params.pos_fnames = pos_fnames(~strcmp(pos_fnames,'NA'));
        %main function for data processing
        [tracks,times] = time_series_analysis_syc(analysis_params, usedLED(ph));
        display('ready to assign')
        all_tracks.(phase) = tracks;
        all_times.(phase) = times + shift_timing(ph)
   end
   all_tracks_vec{jj} = all_tracks; % problem here
   all_times_vec{jj} = all_times;
end
%store data
save([fname_save],'all_times_vec','all_tracks_vec','posvec')

return