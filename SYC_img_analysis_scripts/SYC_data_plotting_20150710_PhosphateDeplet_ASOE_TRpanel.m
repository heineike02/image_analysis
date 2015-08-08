function all_tracks = data_plotting_template()
%
clear
profile off
ipdir = 'C:\Users\susanychen\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

% where the processed data is stored
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\' 
%input information from experiment here
species = 'SC' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,10]    

%std thresh for calling a peak to find a cell
%0.16 seemed to work best for k.lactis on 16JUL images
std_threshSP.KL = 0.16;
std_threshSP.SC = 0.2;

%Tracking parameters
%estimated maximum number of tracks
maxdisp_1x = 4;

% optical amplification 
% 1x or 1.5x
op_amp =  '1x'; %'1.5x'
storeim = 1;


%Strain BMH42 - SC.MSN2-RFP,  KL.MSN2-YFP 
%wellvecSP.SC = {'A9','B9','C9','D9','E9','F9','G9','H9'};


%% Filter single cells traces
% at least 50% of max length

% load the mat files
fname_save = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

% define variables
phases = {'Pre','Post'};
field = 3; % 3 = 'nf'
Xper = 0.95;

% loop through each well
for ii = 1:length(all_tracks_vec)
    all_tracks = all_tracks_vec{ii};
    all_times = all_times_vec{ii};
    %ii 
    % loop through each phase
    for jj = 1:length(phases)
        
        clear tracks_cell; clear desired_track_cell;
        %jj
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
        % observe histograms
        %if jj == 2
        %    figure; histogram(track_lengths(indc_of_tracks));
        %end   
    end
    
    all_tracks_filt_vec{ii} = all_tracks_filt;
    all_times_filt_vec{ii} = all_times_filt;
end
        
% save these filtered traces
fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filtered.mat';
save([fname_save],'all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Mean Time Traces - SYC
figure(1)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

fname_save = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

phases = {'Pre', 'Post'};
channel = ''
color_val = 'g'; %cmap(jj,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
val_to_plot = 'nf';
Nwells = length(all_tracks_vec)

num_subplots = ceil(sqrt(Nwells));

TR_names = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
    'SYC1-70 crz1', 'SYC1-62 Dot6', ...
    'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
    'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
    'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
    'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:length(all_tracks_vec) % loop through each well
    
    subplot(num_subplots,num_subplots,theWells);
            
    all_tracks = all_tracks_vec{theWells};
    all_times = all_times_vec{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(TR_names{theWells})
        axis([0 130 1 4])
    end
    
end

% THIS DOESN'T WORK
% combine pre and post
time_val_cell_pre = timevals_cell(:,1);
time_val_cell_post = timevals_cell(:,2);
time_val_cell = cellfun(@vertcat,time_val_cell_pre,time_val_cell_post, 'UniformOutput', false);
mean_val_cell_pre = mean_val_cell(:,1);
mean_val_cell_post = mean_val_cell(:,2);
mean_val_cell = cellfun(@vertcat,mean_val_cell_pre,mean_val_cell_post, 'UniformOutput', false);
std_val_cell_pre = std_val_cell(:,1);
std_val_cell_post = std_val_cell(:,2);
std_val_cell = cellfun(@vertcat,std_val_cell_pre,std_val_cell_post, 'UniformOutput', false);

%% Mean Time Traces - filtered SYC
figure(2)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

fname_save = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filtered.mat';
load([base_dir,fname_save],'all_times_filt_vec','all_tracks_filt_vec','posvec')

phases = {'Pre', 'Post'};
channel = ''
color_val = 'g'; %cmap(jj,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
val_to_plot = 'nf';
Nwells = length(all_tracks_filt_vec)

num_subplots = ceil(sqrt(Nwells));

TR_names = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
    'SYC1-70 crz1', 'SYC1-62 Dot6', ...
    'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
    'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
    'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
    'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:length(all_tracks_filt_vec) % loop through each well
    
    subplot(num_subplots,num_subplots,theWells);
            
    all_tracks = all_tracks_filt_vec{theWells};
    all_times = all_times_filt_vec{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(TR_names{theWells})
        axis([0 130 1 4])
    end
    
end

% THIS DOESN'T WORK
% combine pre and post
time_val_cell_pre = timevals_cell(:,1);
time_val_cell_post = timevals_cell(:,2);
time_val_cell = cellfun(@vertcat,time_val_cell_pre,time_val_cell_post, 'UniformOutput', false);
mean_val_cell_pre = mean_val_cell(:,1);
mean_val_cell_post = mean_val_cell(:,2);
mean_val_cell = cellfun(@vertcat,mean_val_cell_pre,mean_val_cell_post, 'UniformOutput', false);
std_val_cell_pre = std_val_cell(:,1);
std_val_cell_post = std_val_cell(:,2);
std_val_cell = cellfun(@vertcat,std_val_cell_pre,std_val_cell_post, 'UniformOutput', false);

%% SC.MSN2-RFP plot individual cells 
%% Note the variance of single cell traces
% - variance of the time trace (captures presence of responders and non and
% different types of profiles)
figure(2)
clf
clear t_singlecell_cell; clear y_singlecell_cell;
hold on
channel = ''; %'RFP'
Nwells = length(all_tracks_filt_vec)
num_subplots = ceil(sqrt(Nwells));

TR_names = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
    'SYC1-70 crz1', 'SYC1-62 Dot6', ...
    'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
    'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
    'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
    'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};

%already loaded

%perm = [1,2] ;
%N_RFP = length(perm);
%legend_vec_RFP = legend_vec_RFP(perm);
%cmap = cmap_RFP(perm,:);

for jj = 1:Nwells
    jj
    all_tracks = all_tracks_filt_vec{jj};
    all_times = all_times_filt_vec{jj};
    %subplot(1,2,jj)
    subplot(num_subplots,num_subplots,jj);
    for ph = 1: length(phases) % loop through pre and post
        
        tracks = all_tracks.(phases{ph}); size(tracks)
        timevals = all_times.(phases{ph}); %size(timevals)
        [p] = plot_individual_values(timevals,tracks,channel,'nf');
        hold on;
        
        t_singlecell_cell{jj,ph} = timevals;
        y_singlecell_cell{jj,ph} = tracks;
    end
    axis([0,130,0,10])
    title(TR_names{jj})
end

% can't actually combine the Pre and Post for the single cell data
%% Cluster the single cell traces into a few groups
Nwells_singleCell = y_singlecell_cell(:,2);
Nwells_t_singleCell = t_singlecell_cell(:,2);
    singleCell = Nwells_singleCell{2};
    tsingleCell = Nwells_t_singleCell{2};
    
    
    for tt = 1:length(singleCell)
        [Dist(tt),D,k(tt),w]=dtw(singleCell(1).nf',singleCell(tt).nf');
    end
    
    [idx,C]=kmeans(Dist./k,2);
    ind1=find(idx==1);
    ind2=find(idx==2);
    ind3=find(idx==3);
    ind4=find(idx==4);
    
    figure(1); hold on;
    for yy = 1:length(ind1)
        plot(singleCell(ind1(yy)).nf)
    end
    axis([0 60 1 7])
    
    figure(2); hold on;
    for yy = 1:length(ind2)
        plot(singleCell(ind2(yy)).nf)
    end
    axis([0 60 1 7])

    figure(3); hold on;
    for yy = 1:length(ind3)
        plot(singleCell(ind3(yy)).nf)
    end
    axis([0 60 1 7])
    
    figure(4); hold on;
    for yy = 1:length(ind4)
        plot(singleCell(ind4(yy)).nf)
    end
    axis([0 60 1 4])
    
%% Function that finds the cell associated with trace and makes a movie of that cell
img_dir = '\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20150710_PhosphateDeplet_ASOE_TRpanel\phosphateDeplet\';
Nwells_singleCell = y_singlecell_cell(:,2);
Nwells_t_singleCell = t_singlecell_cell(:,2);

for ii = 3:3%length(Nwells_singleCell{ii}) % loop through each well
    
    singleCell = Nwells_singleCell{ii};
    tsingleCell = Nwells_t_singleCell{ii};
    for kk = 3:3%length(singleCells)
        
  
        img = imread(strcat(img_dir,'RFP_p7_t1.tiff'));
        cX = singleCell(10).Cxloc;
        cY = singleCell(10).Cyloc;
        figure; imshow(img,[]); hold on; plot(cY(1),cX(1),'r.');
        %for jj = 1:length(singleCell(kk).Cxloc)
        %    img = imread(strcat(img_dir,filenames(jj).name));
        %    cX = singleCell(kk).Cxloc(jj);
        %    cY = singleCell(kk).Cyloc(jj);
        %    figure(jj); imshow(img,[]); hold on; plot(cY,cX,'r.');
        %end
    end
end
%% Feature Extraction and Smoothing (Optional)
Nwells_singleCell = y_singlecell_cell(:,2);
Nwells_t_singleCell = t_singlecell_cell(:,2);

var_to_plot = 'nf';
smooth_flag = 1;

% loop through each well
for ii = 1:length(Nwells_singleCell)
    singleCells = Nwells_singleCell{ii};
    tsingleCells = Nwells_t_singleCell{ii};
    
    length(singleCells)
    ii
    % loop through each single cell trace
    for kk = 1:length(singleCells)  
        %kk
        singleCells1 = singleCells(kk);
        % extract features: number of peaks, max height(s), time(s) of max height, rise time (on slope), duration of pulse
        out_Features(kk).cell = extract_singlecell_features_Analysis(tsingleCells, singleCells1,var_to_plot, smooth_flag);
    end
    out_Features_cell(ii).out_Features = out_Features;
%     out_Features_cell(ii).var_over_time =1; 
%     out_Features_cell(ii).numPeaksVar =1;
%     out_Features_cell(ii).maxHeightVar = 1;
%     out_Features_cell(ii).timeMaxHeightVar = 1;
%     out_Features_cell(ii).onSlope
    clear out_Features;
end

%% Plotting histogram distributions of features
%% Plotting correlations among features

%%
% the labels that go with each well
TR_names = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
    'SYC1-70 crz1', 'SYC1-62 Dot6', ...
    'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
    'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
    'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
    'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};
stress_type = {'20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'
    '20150710pHbasic','20150710pHbasic'
    '20150710pHbasic','20150710pHbasic'
    '20150710pHbasic','20150710pHbasic'
    '20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%
% save plotting values/features
fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720.mat';
save([fname_save],'TR_names','stress_type',...
    'time_val_cell','mean_val_cell','std_val_cell',...
    't_singlecell_cell', 'y_singlecell_cell', ...
    'out_Features_cell')