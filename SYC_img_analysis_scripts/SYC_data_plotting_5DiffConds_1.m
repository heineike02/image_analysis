% SYC_data_plotting_FiveDifferentConditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OSMO
% This function:
%     1) Filters single cell traces
%     2) Extracts single cell features
%     3) and Plots summary metric for single cell features
%
% this function takes in 1 Experiment and processes the samples in the
% experiment
%
% %%%% = places where you need to change the script

% set paths
clear
profile off
%%%% ipdir = where your github code is located %%%%
ipdir = '/Users/susanychen/Documents/image_analysis/';
%ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% phases = subfolders in the same experiment %%%%
phases =  {'Pre','Post'}; %,'Post_p2'}
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = '20150709_SorbitolStress_ASOE_TRpanel_20150720.mat';

%% Filter out single cell traces that are too short

% load the mat files
load([base_dir,fname_load],'all_times_vec','all_tracks_vec','posvec')

%%%% field => 'nf' = 3 or 'nmi' = 4 %%%%
field = 3;
%%%% Xper = percent of max trace length %%%%
Xper = 0.95;

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
        
%%%% ffolder_save = folder where filtered traces will be stored %%%%
ffolder_save = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% fname_save = file name of filtered traces %%%%
fname_save = '20150709_SorbitolStress_ASOE_TRpanel_20150720_filter.mat';
save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% folder where the filtered traces live %%%%
ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
fname_load = '20150709_SorbitolStress_ASOE_TRpanel_20150720_filter.mat';
load(strcat(ffolder_load, fname_load),'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
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
        peakFindr.Sel = 0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
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

%%%% TR_names_long =  the labels that go with each well %%%%
TR_names_long = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
    'SYC1-70 crz1', 'SYC1-62 Dot6', ...
    'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
    'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
    'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
    'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};
%%%% TR_names_short =  the labels that go with each well %%%%
TR_names_short = {'Cad1', 'Com2', ...
    'crz1', 'Dot6', ...
    'Maf1', 'Msn2', ...
    'Yap1', 'Stb3', ...
    'Sko1', 'Rtg3', ...
    'Pho4', 'Nrg2','Msn4'};
%%%% stress_type = labels of stress type %%%%
stress_type = {'20150710 sorbStress'}; %,'20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = '/Users/susanychen/Downloads/TempMatlabData20150808/20150709_SorbitolStress_ASOE_TRpanel_20150720_filter.mat';
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names_long', 'TR_names_short','stress_type')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GD2
% This function:
%     1) Filters single cell traces
%     2) Extracts single cell features
%     3) and Plots summary metric for single cell features
%
% this function takes in 1 Experiment and processes the samples in the
% experiment
%
% %%%% = places where you need to change the script

% set paths
clear
profile off
%%%% ipdir = where your github code is located %%%%
ipdir = '/Users/susanychen/Documents/image_analysis/';
%ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% phases = subfolders in the same experiment %%%%
phases =  {'Pre','Post'}; %,'Post_p2'}
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = '20150629_GD_doseresp_Msn2Msn4Maf1Stb3_20150720.mat';

%% Filter out single cell traces that are too short

% load the mat files
load([base_dir,fname_load],'all_times_vec','all_tracks_vec','posvec')

%%%% field => 'nf' = 3 or 'nmi' = 4 %%%%
field = 3;
%%%% Xper = percent of max trace length %%%%
Xper = 0.95;

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
        
%%%% ffolder_save = folder where filtered traces will be stored %%%%
ffolder_save = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% fname_save = file name of filtered traces %%%%
fname_save = '20150629_GD_doseresp_Msn2Msn4Maf1Stb3_20150720_filter.mat';
save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% folder where the filtered traces live %%%%
ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
fname_load = '20150629_GD_doseresp_Msn2Msn4Maf1Stb3_20150720_filter.mat';
load(strcat(ffolder_load, fname_load),'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
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
        peakFindr.Sel = 0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
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

%%%% TR_names_long =  the labels that go with each well %%%%
TR_names_long = {'SYC1-72 Msn2 1%','SYC1-74 Msn4 1%', 'SYC1-75 Maf1 1%', 'SYC1-76 Stb3 1%',...
    'SYC1-72 Msn2 0.5%','SYC1-74 Msn4 0.5%', 'SYC1-75 Maf1 0.5%', 'SYC1-76 Stb3 0.5%',...
    'SYC1-72 Msn2 0.05%','SYC1-74 Msn4 0.05%', 'SYC1-75 Maf1 0.05%', 'SYC1-76 Stb3 0.05%'};
%%%% TR_names_short =  the labels that go with each well %%%%
TR_names_short = {'Msn2 1%','Msn4 1%', 'Maf1 1%', 'Stb3 1%',...
    'Msn2 0.5%','Msn4 0.5%', 'Maf1 0.5%', 'Stb3 0.5%',...
    'Msn2 0.05%','Msn4 0.05%', 'Maf1 0.05%', 'Stb3 0.05%'};
%%%% stress_type = labels of stress type %%%%
stress_type = {'20150629 Glucose Dropout'}; %,'20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = '/Users/susanychen/Downloads/TempMatlabData20150808/20150629_GD_doseresp_Msn2Msn4Maf1Stb3_20150720_filter.mat';
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names_long', 'TR_names_short','stress_type')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GD1
% This function:
%     1) Filters single cell traces
%     2) Extracts single cell features
%     3) and Plots summary metric for single cell features
%
% this function takes in 1 Experiment and processes the samples in the
% experiment
%
% %%%% = places where you need to change the script

% set paths
clear
profile off
%%%% ipdir = where your github code is located %%%%
ipdir = '/Users/susanychen/Documents/image_analysis/';
%ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% phases = subfolders in the same experiment %%%%
phases =  {'Pre','Post'}; %,'Post_p2'}
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = '20150630_GD_binary_BMH35adhOE_20150720.mat';

%% Filter out single cell traces that are too short

% load the mat files
load([base_dir,fname_load],'all_times_vec','all_tracks_vec','posvec')

%%%% field => 'nf' = 3 or 'nmi' = 4 %%%%
field = 3;
%%%% Xper = percent of max trace length %%%%
Xper = 0.95;

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
        
%%%% ffolder_save = folder where filtered traces will be stored %%%%
ffolder_save = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% fname_save = file name of filtered traces %%%%
fname_save = '20150630_GD_binary_BMH35adhOE_20150720_filter.mat';
save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% folder where the filtered traces live %%%%
ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
fname_load = '20150630_GD_binary_BMH35adhOE_20150720_filter.mat';
load(strcat(ffolder_load, fname_load),'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
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
        peakFindr.Sel = 0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
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

%%%% TR_names_long =  the labels that go with each well %%%%
TR_names_long = {'SYC1-62 Dot6', 'SYC1-71 Pho4', 'SYC1-72 Msn2', 'SYC1-65 Rtg3', ...
    'SYC1-67 Cad1' ,'SYC1-68 Nrg2', ...
    'SYC1-69 Sko1', 'SYC1-73 Com2', 'SYC1-64 Yap1'};
%%%% TR_names_short =  the labels that go with each well %%%%
TR_names_short = {'Dot6', 'Pho4', 'Msn2', 'Rtg3', ...
    'Cad1' ,'Nrg2', ...
    'Sko1', 'Com2', 'Yap1'};
%%%% stress_type = labels of stress type %%%%
stress_type = {'20150630 Glucose Dropout'}; %,'20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = '/Users/susanychen/Downloads/TempMatlabData20150808/20150630_GD_binary_BMH35adhOE_20150720_filter.mat';
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names_long', 'TR_names_short','stress_type')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcofluor Zymolyase
% This function:
%     1) Filters single cell traces
%     2) Extracts single cell features
%     3) and Plots summary metric for single cell features
%
% this function takes in 1 Experiment and processes the samples in the
% experiment
%
% %%%% = places where you need to change the script

% set paths
clear
profile off
%%%% ipdir = where your github code is located %%%%
ipdir = '/Users/susanychen/Documents/image_analysis/';
%ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% phases = subfolders in the same experiment %%%%
phases =  {'Pre','Post'}; %,'Post_p2'}
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = '20150716_Calcofluor_Zymolyase_ASOE_5TRs.mat';

%% Filter out single cell traces that are too short

% load the mat files
load([base_dir,fname_load],'all_times_vec','all_tracks_vec','posvec')

%%%% field => 'nf' = 3 or 'nmi' = 4 %%%%
field = 3;
%%%% Xper = percent of max trace length %%%%
Xper = 0.95;

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
        
%%%% ffolder_save = folder where filtered traces will be stored %%%%
ffolder_save = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% fname_save = file name of filtered traces %%%%
fname_save = '20150716_Calcofluor_Zymolyase_ASOE_5TRs_filter.mat';
save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% folder where the filtered traces live %%%%
ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
fname_load = '20150716_Calcofluor_Zymolyase_ASOE_5TRs_filter.mat';
load(strcat(ffolder_load, fname_load),'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
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
        peakFindr.Sel = 0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
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

%%%% TR_names_long =  the labels that go with each well %%%%
TR_names_long = {'SYC1-62 Dot6 calco', 'SYC1-71 Pho4 calco', 'SYC1-72 Msn2 calco', 'SYC1-74 Msn4 calco', 'SYC1-75 Maf1 calco', ...
    'SYC1-75 Maf1 zymo', 'SYC1-74 Msn4 zymo', 'SYC1-72 Msn2 zymo', 'SYC1-71 Pho4 zymo', 'SYC1-62 Dot6 zymo'};
%%%% TR_names_short =  the labels that go with each well %%%%
TR_names_short = {'Dot6 calco', 'Pho4 calco', 'Msn2 calco', 'Msn4 calco', 'Maf1 calco', ...
    'Maf1 zymo', 'Msn4 zymo', 'Msn2 zymo', 'Pho4 zymo', 'Dot6 zymo'};
%%%% stress_type = labels of stress type %%%%
stress_type = {'20150716 Calcofluor Zymolyase'}; %,'20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = '/Users/susanychen/Downloads/TempMatlabData20150808/20150716_Calcofluor_Zymolyase_ASOE_5TRs_filter.mat';
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names_long', 'TR_names_short','stress_type')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pH
% This function:
%     1) Filters single cell traces
%     2) Extracts single cell features
%     3) and Plots summary metric for single cell features
%
% this function takes in 1 Experiment and processes the samples in the
% experiment
%
% %%%% = places where you need to change the script

% set paths
clear
profile off
%%%% ipdir = where your github code is located %%%%
ipdir = '/Users/susanychen/Documents/image_analysis/';
%ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% phases = subfolders in the same experiment %%%%
phases =  {'Pre','Post'}; %,'Post_p2'}
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = '20150716_PhosphateDeplet_pH_ASOE_5TRs_20150720.mat';

%% Filter out single cell traces that are too short

% load the mat files
load([base_dir,fname_load],'all_times_vec','all_tracks_vec','posvec')

%%%% field => 'nf' = 3 or 'nmi' = 4 %%%%
field = 3;
%%%% Xper = percent of max trace length %%%%
Xper = 0.95;

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
        
%%%% ffolder_save = folder where filtered traces will be stored %%%%
ffolder_save = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% fname_save = file name of filtered traces %%%%
fname_save = '20150716_PhosphateDeplet_pH_ASOE_5TRs_20150720_filter.mat';
save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% folder where the filtered traces live %%%%
ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
fname_load = '20150716_PhosphateDeplet_pH_ASOE_5TRs_20150720_filter.mat';
load(strcat(ffolder_load, fname_load),'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
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
        peakFindr.Sel = 0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
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

%%%% TR_names_long =  the labels that go with each well %%%%
TR_names_long = {'SYC1-62 Dot6 phospho', 'SYC1-71 Pho4 phospho', 'SYC1-72 Msn2 phospho', 'SYC1-74 Msn4 phospho', 'SYC1-75 Maf1 phospho', ...
    'SYC1-75 Maf1 pH5.15', 'SYC1-74 Msn4 pH5.15', 'SYC1-72 Msn2 pH5.15', 'SYC1-71 Pho4 pH5.15', 'SYC1-62 Dot6 pH5.15'};
%%%% TR_names_short =  the labels that go with each well %%%%
TR_names_short = {'Dot6 phospho', 'Pho4 phospho', 'Msn2 phospho', 'Msn4 phospho', 'Maf1 phospho', ...
    'Maf1 pH5.15', 'Msn4 pH5.15', 'Msn2 pH5.15', 'Pho4 pH5.15', 'Dot6 pH5.15'};
%%%% stress_type = labels of stress type %%%%
stress_type = {'20150716 pH'}; %,'20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = '/Users/susanychen/Downloads/TempMatlabData20150808/20150716_PhosphateDeplet_pH_ASOE_5TRs_20150720_filter.mat';
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names_long', 'TR_names_short','stress_type')

%% ND

% This function:
%     1) Filters single cell traces
%     2) Extracts single cell features
%     3) and Plots summary metric for single cell features
%
% this function takes in 1 Experiment and processes the samples in the
% experiment
%
% %%%% = places where you need to change the script

% set paths
clear
profile off
%%%% ipdir = where your github code is located %%%%
ipdir = '/Users/susanychen/Documents/image_analysis/';
%ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% phases = subfolders in the same experiment %%%%
phases =  {'Pre','Post'}; %,'Post_p2'}
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = '20150717_AmmoniumSulfateStarvation_ASOE_5TRs_20150720.mat';

%% Filter out single cell traces that are too short

% load the mat files
load([base_dir,fname_load],'all_times_vec','all_tracks_vec','posvec')

%%%% field => 'nf' = 3 or 'nmi' = 4 %%%%
field = 3;
%%%% Xper = percent of max trace length %%%%
Xper = 0.95;

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
        
%%%% ffolder_save = folder where filtered traces will be stored %%%%
ffolder_save = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% fname_save = file name of filtered traces %%%%
fname_save = '20150717_AmmoniumSulfateStarvation_ASOE_5TRs_20150720_filter.mat';
save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% folder where the filtered traces live %%%%
ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
fname_load = '20150717_AmmoniumSulfateStarvation_ASOE_5TRs_20150720_filter.mat';
load(strcat(ffolder_load, fname_load),'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
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
        peakFindr.Sel = 0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
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

%%%% TR_names_long =  the labels that go with each well %%%%
TR_names_long = {'SYC1-62 Dot6 noCSM', 'SYC1-71 Pho4 noCSM', 'SYC1-72 Msn2 noCSM', 'SYC1-74 Msn4 noCSM', 'SYC1-75 Maf1 noCSM', ...
    'SYC1-75 Maf1 CSM', 'SYC1-74 Msn4 CSM', 'SYC1-72 Msn2 CSM', 'SYC1-71 Pho4 CSM', 'SYC1-62 Dot6 CSM'};
%%%% TR_names_short =  the labels that go with each well %%%%
TR_names_short = {'Dot6 noCSM', 'Pho4 noCSM', 'Msn2 noCSM', 'Msn4 noCSM', 'Maf1 noCSM', ...
    'Maf1 CSM', 'Msn4 CSM', 'Msn2 CSM', 'Pho4 CSM', 'Dot6 CSM'};
%%%% stress_type = labels of stress type %%%%
stress_type = {'20150717 Ammonium sulfate depletion'}; %,'20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic'...
    %'20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = '/Users/susanychen/Downloads/TempMatlabData20150808/20150717_AmmoniumSulfateStarvation_ASOE_5TRs_20150720_filter.mat';
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names_long', 'TR_names_short','stress_type')