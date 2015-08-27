function all_tracks = data_plotting_template()
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
fname_load = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720.mat';

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
fname_save = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save(strcat(ffolder_save, fname_save),'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% folder where the filtered traces live %%%%
ffolder_load = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%ffolder_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% filename of filtered traces %%%%
fname_load = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
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
stress_type = {'20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%% -- stopped here
fname_save = '/Users/susanychen/Downloads/TempMatlabData20150808/20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
%fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names_long', 'TR_names_short','stress_type')

%% Mean Time Traces - filtered SYC (filtered, not smoothed)
figure(1)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% specify the .mat file to load %%%%
fname_load = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
%fname_load = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
load([base_dir,fname_load],'all_times_filt_vec','all_tracks_filt_vec','out_Features_cell','TR_names_short', 'stress_type')
%%%% specify the subfolders present in the experiment %%%%
phases = {'Pre', 'Post'};
%%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
channel = '';
%%%% specify color of lines in plot %%%%
color_val = 'g'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};
%%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
val_to_plot = 'nf';

Nwells = length(all_tracks_filt_vec);

num_subplots = ceil(sqrt(Nwells));

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:length(all_tracks_filt_vec) % loop through each well
    
    subplot(num_subplots,num_subplots,theWells);
            
    all_tracks = all_tracks_filt_vec{theWells};
    all_times = all_times_filt_vec{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues_syc(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(TR_names_short{theWells})
        %%%% specify the x and y range of the plot %%%%
        axis([0 130 1 6])
    end
    
end

print('/Users/susanychen/Downloads/TempMatlabData20150808/MeanCellPlot_20150710PhosphateDeplet.eps', '-depsc')
%print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MeanCellPlot.pdf','-dpdf')

%% SC.MSN2-RFP plot individual cells (filtered, not smoothed)
% - variance of the time trace (captures presence of responders and non and
% different types of profiles)
figure(2)
clf
clear t_singlecell_cell; clear y_singlecell_cell;
hold on
%%%% specify channel name %%%%
channel = ''; %'RFP'

Nwells = length(all_tracks_filt_vec);
num_subplots = ceil(sqrt(Nwells));

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
    %%%% specify x and y range of plot %%%%
    axis([0,130,0,10])
    title(TR_names_short{jj})
end

print('/Users/susanychen/Downloads/TempMatlabData20150808/SingleCellPlot_20150710PhosphateDeplet.eps', '-depsc')
%print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\SingleCellPlot.pdf','-dpdf')

%% Plotting "Number of Peaks" Feature
TR_names1 = {'cad1', 'com2', 'crz1', 'dot6', 'maf1', 'msn2' ,'yap1', 'stb3', 'sko1', 'rtg3', 'pho4' ,'nrg2', 'msn4'};

numWells = length(out_Features_cell);
splot = ceil(sqrt(numWells));

clear uniqueVals
clear peakNumVecCell

for jj = 1:numWells
    
    display(strcat('This is jj: ', num2str(jj)))
    
    numCells = length(out_Features_cell(jj).out_Features);
    
    for ii = 1:numCells
        numPeaksVec(ii) = out_Features_cell(jj).out_Features(ii).cell.NumPeaks;
    end
    
    % find unique number of peaks
    [vals, indc] = unique(numPeaksVec);
    
    uniqueVals{jj} = vals;
    peakNumVecCell{jj} = numPeaksVec;

    clear numPeaksVec
end

numDiffPeaks = length(unique(cell2mat(uniqueVals)));
diffPeaks = unique(cell2mat(uniqueVals));

initMatrix = nan(numDiffPeaks, numWells);
for gg = 1:numDiffPeaks
    display(strcat('gg: ',num2str(gg)))
    for hh = 1:numWells
        display(num2str(hh))
        if find(peakNumVecCell{hh} == diffPeaks(gg)) 
            display('exist')
            totCellNum = numel(peakNumVecCell{hh});
            currPeakNum = numel(find(peakNumVecCell{hh} == diffPeaks(gg)));
            propPeakNum = currPeakNum/totCellNum;
            initMatrix(gg,hh) = propPeakNum;
            
        else
            display('dne')
            initMatrix(gg,hh) = NaN;
        end
    end
end

% plot a stacked bar graph
figure(1)
bar(1:numWells, initMatrix', 0.5,'stack') %[initMatrix(1,:)', initMatrix(2,:)', initMatrix(3,:)', initMatrix(4,:)'], 0.5, 'stack')
set(gca, 'xticklabel', TR_names1)
xlabel('TRs'); ylabel('Number of Cells (count)')
legend(strsplit(num2str(diffPeaks)))
title('Proportion of Number of Peaks for TRs')

print('/Users/susanychen/Downloads/TempMatlabData20150808/NumPeaksPlot_20150710PhosphateDeplet.eps', '-depsc')

%% FIRST PEAK ONLY:
% SUBPLOT: MAX HEIGHT, TIME MAX HEIGHT, PULSE WIDTH, OFFSLOPE, ONSLOPE

TR_names1 = {'crz1', 'dot6', 'maf1', 'msn2' ,'yap1', 'stb3', 'sko1', 'rtg3', 'pho4' ,'nrg2', 'msn4'}; 

numWells = length(out_Features_cell);

pp =1;
first = 1;
clear maxHeightVec
clear maxHeightVecCell
clear maxHeightTimeVec
clear maxHeightTimeVecCell
clear pulseWidthVec
clear pulseWidthVecCell
clear offSlopeVec
clear offSlopeVecCell
clear onSlopeVec
clear onSlopeVecCell

for jj = 3:numWells

    numCells = length(out_Features_cell(jj).out_Features);
    for ii = 1:numCells

        % filter out the no peak data
        if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
            
            % only find the first peak
            if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == 1 & first == 1

                maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1);
                maxHeightTimeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1);
                pulseWidthVec(pp) = out_Features_cell(jj).out_Features(ii).cell.PulseWidth(1);
                offSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OffSlope(1);
                onSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OnSlope(1);
                
                %if ii >30 & ii <32;
                %    figure;
                %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
                %    hold on;
                %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
                %end
                
                pp=pp+1;
            end
        end

    end
    
    maxHeightVecCell{jj-2} = maxHeightVec;
    maxHeightTimeVecCell{jj-2} = maxHeightTimeVec;
    pulseWidthVecCell{jj-2} = pulseWidthVec;
    offSlopeVecCell{jj-2} = offSlopeVec;
    onSlopeVecCell{jj-2} = onSlopeVec;
end


% % %             % look at all the peaks
% % %             for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% % %                 maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(kk);
% % %                 pp=pp+1;
% % %             end
            
%print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight.eps','-depsc')

figure(2); 
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',TR_names1, 'yLabel', 'Maximum Height (nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Maximum Height Time (min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Pulse Width (min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Off Slope (nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'On Slope (nuc/cyt per 2min)'); title('On Slope');
suptitle('First Peak Features');

maxHeightVecCell_1 = maxHeightVecCell;
maxHeightTimeVecCell_1 = maxHeightTimeVecCell;
pulseWidthVecCell_1 = pulseWidthVecCell;
offSlopeVecCell_1 = offSlopeVecCell;
onSlopeVecCell_1 = onSlopeVecCell;

print('/Users/susanychen/Downloads/TempMatlabData20150808/FirstPeakFeaturesPlot_20150710PhosphateDeplet.eps', '-depsc')

%% SECOND PEAK ONLY:
% SUBPLOT: MAX HEIGHT, TIME MAX HEIGHT, PULSE WIDTH, OFFSLOPE, ONSLOPE
numWells = length(out_Features_cell);
pp =1;
clear maxHeightVec
clear maxHeightVecCell
clear maxHeightTimeVec
clear maxHeightTimeVecCell
clear pulseWidthVec
clear pulseWidthVecCell
clear offSlopeVec
clear offSlopeVecCell
clear onSlopeVec
clear onSlopeVecCell

for jj = 3:numWells
    numCells = length(out_Features_cell(jj).out_Features);
    for ii = 1:numCells      
            % only find the first peak
            if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == 2 %& first == 2

                maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(2);
                maxHeightTimeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(2);
                pulseWidthVec(pp) = out_Features_cell(jj).out_Features(ii).cell.PulseWidth(2);
                offSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OffSlope(2);
                onSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OnSlope(2);
                
                %if ii >30 & ii <32;
                %    figure;
                %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
                %    hold on;
                %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
                %end
                
                pp=pp+1;
            end
    end
    
    maxHeightVecCell{jj-2} = maxHeightVec;
    maxHeightTimeVecCell{jj-2} = maxHeightTimeVec;
    pulseWidthVecCell{jj-2} = pulseWidthVec;
    offSlopeVecCell{jj-2} = offSlopeVec;
    onSlopeVecCell{jj-2} = onSlopeVec;
end

% % %             % look at all the peaks
% % %             for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% % %                 maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(kk);
% % %                 pp=pp+1;
% % %             end

figure(3); 
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',TR_names1, 'yLabel', 'Maximum Height (nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Maximum Height Time (min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Pulse Width (min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Off Slope (nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'On Slope (nuc/cyt per 2min)'); title('On Slope');
suptitle('Second Peak Features');

maxHeightVecCell_2 = maxHeightVecCell;
maxHeightTimeVecCell_2 = maxHeightTimeVecCell;
pulseWidthVecCell_2 = pulseWidthVecCell;
offSlopeVecCell_2 = offSlopeVecCell;
onSlopeVecCell_2 = onSlopeVecCell;

print('/Users/susanychen/Downloads/TempMatlabData20150808/SecondPeakFeaturesPlot_20150710PhosphateDeplet.eps', '-depsc')

%% COMPARING FIRST AND SECOND PEAK FEATURES - plot over each other
figure(3); 
subplot(2,3,1); 
distributionPlot(maxHeightVecCell_1,'distWidth',0.9, 'color', 'g', 'addSpread',0, 'showMM',2,'xNames',TR_names1, 'yLabel', 'Maximum Height (nuc/cyt)'); title('Maximum Height'); hold on;
distributionPlot(maxHeightVecCell_2,'distWidth',0.9, 'color', 'c', 'addSpread',0, 'showMM',2); legend('peak1', 'peak2')
subplot(2,3,2); 
distributionPlot(maxHeightTimeVecCell_1,'distWidth',0.9,'color', 'g', 'addSpread',0,'showMM',2,'xNames',TR_names1, 'yLabel', 'Maximum Height Time (min)'); title('Maximum Height Time');hold on;
distributionPlot(maxHeightTimeVecCell_2,'distWidth',0.9, 'color', 'c', 'addSpread',0, 'showMM',2);
subplot(2,3,3); 
distributionPlot(pulseWidthVecCell_1,'distWidth',0.9,'color', 'g', 'addSpread',0,'showMM',2,'xNames',TR_names1, 'yLabel', 'Pulse Width (min)'); title('Pulse Width');hold on;
distributionPlot(pulseWidthVecCell_2,'distWidth',0.9, 'color', 'c', 'addSpread',0, 'showMM',2);
subplot(2,3,4); 
distributionPlot(offSlopeVecCell_1,'distWidth',0.9,'color', 'g', 'addSpread',0,'showMM',2,'xNames',TR_names1, 'yLabel', 'Off Slope (nuc/cyt per 2min)'); title('Off Slope');hold on;
distributionPlot(offSlopeVecCell_2,'distWidth',0.9, 'color', 'c', 'addSpread',0, 'showMM',2);
subplot(2,3,5); 
distributionPlot(onSlopeVecCell_1,'distWidth',0.9,'color', 'g', 'addSpread',0,'showMM',2,'xNames',TR_names1, 'yLabel', 'On Slope (nuc/cyt per 2min)'); title('On Slope');hold on;
distributionPlot(onSlopeVecCell_2,'distWidth',0.9, 'color', 'c', 'addSpread',0, 'showMM',2);
suptitle('Comparison of First and Second Peak Features');

print('/Users/susanychen/Downloads/TempMatlabData20150808/FirstAndSecondPeakFeaturesPlot_20150710PhosphateDeplet.eps', '-depsc')

%% Plotting correlations and correlation coefficients

maxHeightVecCell_2 = maxHeightVecCell;
maxHeightTimeVecCell_2 = maxHeightTimeVecCell;
pulseWidthVecCell_2 = pulseWidthVecCell;
offSlopeVecCell_2 = offSlopeVecCell;
onSlopeVecCell_2 = onSlopeVecCell;

% correlations of 1st peak properties
for ee = 1:length(maxHeightVecCell_1)
    
    curr1 = maxHeightVecCell_1{ee};
    curr2 = maxHeightTimeVecCell_1{ee};
    corr_heightVSheightTime1(ee) = corr2(curr1,curr2);
  
    curr2 = pulseWidthVecCell_1{ee};
    corr_heightVSpulseWidth1(ee) = corr2(curr1,curr2);

    curr2 = offSlopeVecCell_1{ee};
    corr_heightVSoffSlope1(ee) = corr2(curr1,curr2);

    curr2 = onSlopeVecCell_1{ee};
    corr_heightVSonSlope1(ee) = corr2(curr1,curr2);
 
end

% correlations of 2nd peak properties
for ee = 1:length(maxHeightVecCell_2)
    
    curr1 = maxHeightVecCell_2{ee};
    curr2 = maxHeightTimeVecCell_2{ee};
    corr_heightVSheightTime2(ee) = corr2(curr1,curr2);
  
    curr2 = pulseWidthVecCell_2{ee};
    corr_heightVSpulseWidth2(ee) = corr2(curr1,curr2);

    curr2 = offSlopeVecCell_2{ee};
    corr_heightVSoffSlope2(ee) = corr2(curr1,curr2);

    curr2 = onSlopeVecCell_2{ee};
    corr_heightVSonSlope2(ee) = corr2(curr1,curr2);

end

figure(1);
subplot(4,1,1)
bar([corr_heightVSheightTime1',corr_heightVSheightTime2']); title('peak height vs peak time'); ylabel('correlation coefficient')
legend('first peak', 'second peak'); xlim([0 12]); set(gca, 'xticklabels', TR_names_short); %rotateXLabels(gca, 90);
subplot(4,1,2);
bar([corr_heightVSpulseWidth1',corr_heightVSpulseWidth2']); title('peak height vs peak width');
xlim([0 12]); set(gca, 'xticklabels', TR_names_short); %rotateXLabels(gca, 90);
subplot(4,1,3);
bar([corr_heightVSoffSlope1',corr_heightVSoffSlope2']); title('peak height vs off slope')
xlim([0 12]); set(gca, 'xticklabels', TR_names_short); %rotateXLabels(gca, 90);
subplot(4,1,4); 
bar([corr_heightVSonSlope1',corr_heightVSonSlope2']); title('peak height vs on slope')
xlim([0 12]); set(gca, 'xticklabels', TR_names_short); %rotateXLabels(gca, 90);
% get the labeling right here %
clear corr_heightVSheightTime1; clear corr_heightVSpulseWidth1; clear corr_heightVSoffSlope1; clear corr_heightVSonSlope1;

clear corr_heightVSheightTime2; clear corr_heightVSpulseWidth2; clear corr_heightVSoffSlope2; clear corr_heightVSonSlope2;

base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\';
save(strcat(base_dir,'OrderOfTRNames.mat'), 'TR_names1')

%% Data Analysis Thoughts %%
% 1. It will be nice to see the traces, peaks, and other demarcations for
% different parts of the distribution
% 2. Run the data plots by Raj
% 3. Okay to look for SIMILAR trends for all TRs for a given input
% 4. But really need to look at how TRs respond in different environmental
% conditions
% 5. consider plotting things WRT Msn2