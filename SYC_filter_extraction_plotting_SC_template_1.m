function [] = SYC_filter_extraction_plotting_SC_template_1(ipdir, base_dir, phases, fname_load, well_order, perturbation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function:
%     1) Filters single cell traces
%     2) Extracts single cell features
%     3) and Plots summary metric for single cell features
%
% this function takes in 1 Experiment and processes the samples in the
% experiment
%
% %%%% = places where you need to change the script
%
% Inputs to function:
%       1. ipdir = directory where github image analysis code is stored
%       2. base_dir = directory where the processed .mat file of single
%       cell traces are stored (from previous processing scripts)
%       3. phases = the "before" and "after" environmental conditions in a
%       given experiment
%       4. fname_load = filename of the processed .mat file of single cell
%       traces
%
%
% Outputs to function:
%       there are no outputs
%
% Example of inputs:
%       1. ipdir = '/Users/susanychen/Documents/image_analysis/';
%       2. base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%       3. phases = {'Pre','Post'};
%       4. fname_load= '20150709_SorbitolStress_ASOE_TRpanel_20150720';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set paths
clear
profile off
ipdir = ipdir;
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = base_dir;
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% phases = subfolders in the same experiment %%%%
phases = phases; 
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = fname_load;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter out single cell traces that are too short

% load the mat files
load([base_dir,fname_load],'all_times_vec','all_tracks_vec','posvec')

%%%% field => 'nf' = 3 or 'nmi' = 4 %%%%
field = 3;
%%%% Xper = percent of max trace length %%%%
Xper = 0.9;
%Xper = 0.95;

% loop through each well
for ii = 1:length(all_tracks_vec)
    all_tracks = all_tracks_vec{ii};
    all_times = all_times_vec{ii};

    % loop through each phase
    for jj = 1:length(phases)
        
        clear tracks_cell; clear desired_track_cell;

        tracks = all_tracks.(phases{jj});
        timevals = all_times.(phases{jj});
        
        % filteres out traces that are < X% of max trace length
        filtered_tracks_cell = filter_singlecell_Analysis(tracks, field, Xper);
        
        % convert back to structure from cell
        filtered_tracks_struct = cell2struct(filtered_tracks_cell, {'Cxloc','Cyloc', 'nf', 'nmi', 'times', 'length', 'pos'}, 1);
        
        % build back the storage structure
        if jj == 1 % this is Pre
            all_tracks_filt.Pre = filtered_tracks_struct;
            all_times_filt.Pre = timevals;
        elseif jj == 2 % this is Post
            all_tracks_filt.Post = filtered_tracks_struct;
            all_times_filt.Post = timevals;
        end
        
        % visual check
        %if jj == 2
        %    figure; histogram(track_lengths(indc_of_tracks));
        %end  
        
    end
    
    all_tracks_filt_vec{ii} = all_tracks_filt;
    all_times_filt_vec{ii} = all_times_filt;
end

cell_num_filt_vec = [];
for hh = 1:length(all_tracks_filt_vec) 
    cell_num_filt_vec(hh) = length(all_tracks_filt_vec{hh}.Post);
end

% How many cells remain after initial filtering?
%%%% ffolder_save = folder where filtered traces will be stored %%%%
%ffolder_save = base_dir;
%ffolder_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% fname_save = file name of filtered traces %%%%
%fname_save = [fname_load,'_filter.mat'];
%save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec', 'cell_num_filt_vec')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature Extraction and (Optional) Smoothing
%%%% folder where the filtered traces live %%%%
%ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
%fname_load_1 = [fname_load, 'filter.mat'];
%load(strcat(ffolder_load, fname_load_1),'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
clear out_Features_cell;
clear t_singlecell_cell; clear y_singlecell_cell;

% set variables
%%%% If only 1 channel, leave blank, if more than 1 channel, specify the
%%%% channel name %%%%
channel = ''; %'RFP'
Nwells = length(all_tracks_filt_vec);

% converts structure of single cell traces from struct to cell array
for jj = 1:Nwells
    all_tracks = all_tracks_filt_vec{jj};
    all_times = all_times_filt_vec{jj};

    for ph = 1: length(phases) % loop through pre and post
        
        tracks = all_tracks.(phases{ph}); 
        timevals = all_times.(phases{ph});
        
        t_singlecell_cell{jj,ph} = timevals;
        y_singlecell_cell{jj,ph} = tracks;
    end
end

Nwells_singleCell = y_singlecell_cell(:,2);
Nwells_t_singleCell = t_singlecell_cell(:,2);

% set the variables for feature extraction
%%%% 'nf' = nuclear localization %%%%
var_to_plot = 'nf';
%%%% 1 = smooth, 0 = don't smooth %%%%
smoothFn.flag = 1;
%%%% method for smoothing, see smooth function for more options %%%%
smoothFn.method = 'sgolay';
%%%% the window of smoothing, value found empirically %%%%
smoothFn.span = 0.4; % span of 40%

% Smoothes and extracts features from single cell traces
% loop through each well
for ii = 1:length(Nwells_singleCell)
    singleCells = Nwells_singleCell{ii};
    tsingleCells = Nwells_t_singleCell{ii};
    
    length(singleCells)
    
    display(strcat('This is well: ', num2str(ii)))
    
    % loop through each single cell trace
    for kk = 1:length(singleCells)  
         
        singleCells1 = singleCells(kk);
        
        display(strcat('This is cell: ', num2str(kk)))
        
        %%%% threshold for peakfinder function, empirically derived %%%%
        peakFindr.Sel = 0.1; %0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
        %%%% in peakfinder function, lower bound of peak %%%%
        peakFindr.Thresh = 1.5;
        %%%% empirically derived, only change if want to lower the threshold %%%%
        peakFindr.StdMinus = 0; %0.2; 
        
        % test plot
        %figure; plot(singleCells1.nf);
        %%plot(out_Features(kk).cell.TimeMaxHeight', out_Features(kk).cell.MaxHeight, 'ro');
        
        % extract features: number of peaks, max height(s), time(s) of max height, rise time (on slope, off slope, on slope time, off slope time), duration of pulse
        out_Features(kk).cell = extract_singlecell_features_Analysis(tsingleCells, singleCells1,var_to_plot, smoothFn, peakFindr);
        
        % test plot
        %figure; plot(out_Features(kk).cell.TimeTraceTimes',out_Features(kk).cell.TimeTrace); hold on;
        %plot(out_Features(kk).cell.TimeMaxHeight', out_Features(kk).cell.MaxHeight, 'ro');

        %if isfield(out_Features(kk).cell, 'OffSlopeTime')
        %    plot(out_Features(kk).cell.OffSlopeTime,2*ones(1,length(out_Features(kk).cell.OffSlopeTime)),'c.');
        %end
        %if isfield(out_Features(kk).cell, 'OnSlopeTime')
        %    plot(out_Features(kk).cell.OnSlopeTime, 2*ones(1,length(out_Features(kk).cell.OffSlopeTime)),'c.')
        %end
        %axis([0 140 1 10])  
         
    end
    
    % store features
    out_Features_cell(ii).out_Features = out_Features;

    clear out_Features;
end

% number of cells with peaks and no peaks
numWells = length(out_Features_cell);
zeroPeakNum = [];
hasPeakNum = [];
for gg = 1:numWells
    numCells = length(out_Features_cell(gg).out_Features);
    ee = 1;
    ff = 1;
for qq = 1:numCells
    currPeakNum = out_Features_cell(gg).out_Features(qq).cell.NumPeaks;
    if currPeakNum == 0;    
        ee = ee+1;
    else  
        ff = ff+1;
    end
end
    out_Features_cell(gg).zeroPeakNum = ee;
    out_Features_cell(gg).hasPeakNum = ff;
end

% %%%% TR_names_long =  the labels that go with each well %%%%
% TR_names_long = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
%     'SYC1-70 crz1', 'SYC1-62 Dot6', ...
%     'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
%     'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
%     'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
%     'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};
% %%%% TR_names_short =  the labels that go with each well %%%%
% TR_names_short = {'Cad1', 'Com2', ...
%     'crz1', 'Dot6', ...
%     'Maf1', 'Msn2', ...
%     'Yap1', 'Stb3', ...
%     'Sko1', 'Rtg3', ...
%     'Pho4', 'Nrg2','Msn4'};
%%%% stress_type = labels of stress type %%%%
%stress_type = {'20150710 sorbStress'};
well_order = well_order;
perturbation = perturbation;

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = [base_dir,fname_load,'_filter.mat'];
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'well_order', 'perturbation', 'cell_num_filt_vec')
