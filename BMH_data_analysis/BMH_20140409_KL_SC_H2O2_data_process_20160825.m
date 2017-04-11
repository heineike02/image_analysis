function all_tracks = BMH_20140409_KL_SC_H2O2_data_process_20160825()
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
fname_conv = 'JSO';
analysis_params.fname_conv = fname_conv;

% ?? doublecheck this
analysis_params.storeim = 1;

%Base directory for image files
%base_dir = '/Users/liamholt/Documents/MATLAB/0-CRZ1-test/20150527-CRZ1-noCa2_w1488 TIRF low_'
base_dir = '\\elsamad.ucsf.edu\Data\Scratch\Ben\20140409\'
%filenames for each phase of imaging.  If only one phase of imaging code
%will still work - this basically is the directory that Micromanager uses
%to store all image files from a certain set of images. 
%example: 
%
%phase = {'Exp'}
%imdirPhase.Exp = [base_dir,'Exp',filesep]
%

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
%These are the absolute times at which each phase starts.
shift_timing = [0,4]    

%imdirPhase.Exp = [base_dir,'Exp',filesep]
imdirPhase.Pre = [base_dir,'TM_Sorb_H2O2_pre',filesep];
imdirPhase.Post = [base_dir,'TM_Sorb_H2O2_post',filesep];

species_list = {'KL','SC'}
%species_circ = {'circKL_1p5x.mat','circSC_1p5x.mat'}
species_fname_save = {'processed_data_20140409_H202_KL','processed_data_20140409_H202_SC' }

%% Plate information
%Pos 1-2: SC Cont 
%Pos 3-4: KL YFP Cont 
%Pos 5-6: KL RFP Cont
%Pos 7-8: SC TM 5ug/ml
%Pos 9-10: KL YFP TM 5ug/ml
%Pos 11-12: KL RFP TM 5ug/ml 
%Pos 13-14: SC SOrb 0.25M
%Pos 15-16: KL YFP Sorb
%Pos 17-18: KL RFP Sorb
%Pos 19-20: SC H2O2 1.0mM
%Pos 21-22: KL YFP H2O2
%Pos 23-24: KL RFP H2O2

%Positions 11 and 12 had large bright spots and tracking broke at t89
%Many of these don't have enough cells. 
%this does not try to calculate TM perturbation. 

% Using JSO's matlab image code, the folder names are different - here is
% an example of the position vector in that case
species_posvec_KL = {'p5','p6';
 'p17','p18';
 'p23','p24'};

species_posvec_SC = {'p1','p2';
 'p13','p14';
 'p19','p20'};

species_posvec = {species_posvec_KL,species_posvec_SC}

for species_ind = 1:2
    species = species_list{species_ind};
    op_amp = '1p5x';
    maxdisp_1x = 4;
    [circ, siz, rad, maxdisp, std_thresh] = species_magnification_params(species, op_amp, ipdir, maxdisp_1x);
    
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
    
    %in FindCells, after convolution, find possible locations of cells with generous threshold
    %actual threshold is this factor times the mean absolute deviation plus the median.  
    analysis_params.thresh = 2;

    %Tracking Parameter: estimated maximum displacement
    analysis_params.maxdisp = maxdisp;

    %Tracking Parameter: 
    analysis_params.track_memory = 2;

    %Tracking Parameter
    analysis_params.min_points_for_traj = 4;

    %std thresh for calling a peak to find a cell
    analysis_params.std_thresh = std_thresh;

    %total number of channels of images taken
    analysis_params.channels = {'BF','RFP','YFP'};

    %Channels we want to process images from.
    %Use a cell. Can use more than one channel  If just one channel just list it as cell with one value i.e. {'RFP'}. 
    channels_to_image = {'RFP'};  
    analysis_params.channels_to_image = channels_to_image;

    %appends this onto base_dir - right now base_dir has a bit of a filename on it
    %as well. 
    fname_save = species_fname_save{species_ind};


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
    time_calc_phase.Pre =  2.5;
    time_calc_phase.Post = 2.5;
    %time_calc_phase.Post = 'metadata_parsed.txt';





    %% Background image

    bgimg = 0; 
    %1 if you have a background image, 0 if you don't. 
    %Gets background image depending on channel to image
    %Collect background images using micromanager 
    if bgimg == 0  %set default to 1
        for jj = 1:length(channels_to_image)
            imbg.(channels_to_image{jj}) = 1;
        end

    else
        %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
        for jj = 1:length(channels_to_image)
            channel = channels_to_image{jj};
            if strcmp(channel,'RFP')
                 imbg_jj = imread([base_dir,'BG',filesep,'img_000000000_RFP_001.tif']);
            elseif strcmp(channel,'YFP')
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


    posvec = species_posvec{species_ind};

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
    Nconditions = size(posvec,1);
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

end
return
