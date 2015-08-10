%%%% direction1: filter out all the NO PEAKS and then plot the mean/median %%%%%%
%%%% direction2: show what 0 peak, 1 peak, 2 peak, 3 peak looks like %%%%%%

%% This Script plots the same summary plots as for each of the experiments, but across conditions for only 5 TRs

base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
experi_name = '20150629_GD_doseresp_Msn2Msn4Maf1Stb3_20150720_filter';
load([base_dir,experi_name])
%TR_names_long(9:12);
%TR_names_short(9:12);
GD_Msn2_time = all_times_filt_vec(9);
GD_Msn4_time = all_times_filt_vec(10);
GD_Maf1_time = all_times_filt_vec(11);
GD_Msn2_track = all_tracks_filt_vec(9);
GD_Msn4_track = all_tracks_filt_vec(10);
GD_Maf1_track = all_tracks_filt_vec(11);
GD_Msn2_feat = out_Features_cell(9).out_Features;
GD_Msn4_feat = out_Features_cell(10).out_Features;
GD_Maf1_feat = out_Features_cell(11).out_Features;

base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
experi_name = '20150630_GD_binary_BMH35adhOE_20150720_filter';
load([base_dir,experi_name])
GD_Dot6_time = all_times_filt_vec(1);
GD_Pho4_time = all_times_filt_vec(2);
GD_Msn2_1_time = all_times_filt_vec(3);
GD_Dot6_track = all_tracks_filt_vec(1);
GD_Pho4_track = all_tracks_filt_vec(2);
GD_Msn2_1_track = all_tracks_filt_vec(3);
GD_Dot6_feat = out_Features_cell(1).out_Features;
GD_Pho4_feat = out_Features_cell(2).out_Features;
GD_Msn2_1_feat = out_Features_cell(3).out_Features;

base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
experi_name = '20150709_SorbitolStress_ASOE_TRpanel_20150720_filter';
load([base_dir,experi_name])
OS_Dot6_time = all_times_filt_vec(4);
OS_Maf1_time = all_times_filt_vec(5);
OS_Msn2_time = all_times_filt_vec(6);
OS_Pho4_time = all_times_filt_vec(11);
OS_Msn4_time = all_times_filt_vec(13);
OS_Dot6_track = all_tracks_filt_vec(4);
OS_Maf1_track = all_tracks_filt_vec(5);
OS_Msn2_track = all_tracks_filt_vec(6);
OS_Pho4_track = all_tracks_filt_vec(11);
OS_Msn4_track = all_tracks_filt_vec(13);
OS_Dot6_feat = out_Features_cell(4).out_Features;
OS_Maf1_feat = out_Features_cell(5).out_Features;
OS_Msn2_feat = out_Features_cell(6).out_Features;
OS_Pho4_feat = out_Features_cell(11).out_Features;
OS_Msn4_feat = out_Features_cell(13).out_Features;

base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
experi_name = '20150716_Calcofluor_Zymolyase_ASOE_5TRs_filter';
load([base_dir,experi_name])
ZM_Maf1_time = all_times_filt_vec(6);
ZM_Msn4_time = all_times_filt_vec(7);
ZM_Msn2_time = all_times_filt_vec(8);
ZM_Pho4_time = all_times_filt_vec(9);
ZM_Dot6_time = all_times_filt_vec(10);
ZM_Maf1_track = all_tracks_filt_vec(6);
ZM_Msn4_track = all_tracks_filt_vec(7);
ZM_Msn2_track = all_tracks_filt_vec(8);
ZM_Pho4_track = all_tracks_filt_vec(9);
ZM_Dot6_track = all_tracks_filt_vec(10);
ZM_Maf1_feat = out_Features_cell(6).out_Features;
ZM_Msn4_feat = out_Features_cell(7).out_Features;
ZM_Msn2_feat = out_Features_cell(8).out_Features;
ZM_Pho4_feat = out_Features_cell(9).out_Features;
ZM_Dot6_feat = out_Features_cell(10).out_Features;

base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
experi_name = '20150716_PhosphateDeplet_pH_ASOE_5TRs_20150720_filter';
load([base_dir,experi_name])
PH_Dot6_time = all_times_filt_vec(1);
PH_Pho4_time = all_times_filt_vec(2);
PH_Msn2_time = all_times_filt_vec(3);
PH_Msn4_time = all_times_filt_vec(4);
PH_Maf1_time = all_times_filt_vec(5);
PH_Dot6_track = all_tracks_filt_vec(1);
PH_Pho4_track = all_tracks_filt_vec(2);
PH_Msn2_track = all_tracks_filt_vec(3);
PH_Msn4_track = all_tracks_filt_vec(4);
PH_Maf1_track = all_tracks_filt_vec(5);
PH_Dot6_feat = out_Features_cell(1).out_Features;
PH_Pho4_feat = out_Features_cell(2).out_Features;
PH_Msn2_feat = out_Features_cell(3).out_Features;
PH_Msn4_feat = out_Features_cell(4).out_Features;
PH_Maf1_feat = out_Features_cell(5).out_Features;

base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
experi_name = '20150717_AmmoniumSulfateStarvation_ASOE_5TRs_20150720_filter';
load([base_dir,experi_name])
ND_Dot6_time = all_times_filt_vec(1); % could change to 1-5
ND_Pho4_time = all_times_filt_vec(2);
%ND_Msn2_time = all_times_filt_vec(8);
ND_Msn4_time = all_times_filt_vec(4);
ND_Maf1_time = all_times_filt_vec(5);
ND_Dot6_track = all_tracks_filt_vec(1);
ND_Pho4_track = all_tracks_filt_vec(2);
%ND_Msn2_track = all_tracks_filt_vec(8);
ND_Msn4_track = all_tracks_filt_vec(4);
ND_Maf1_track = all_tracks_filt_vec(5);
ND_Dot6_feat = out_Features_cell(1).out_Features;
ND_Pho4_feat = out_Features_cell(2).out_Features;
%ND_Msn2_feat = out_Features_cell(8).out_Features;
ND_Msn4_feat = out_Features_cell(4).out_Features;
ND_Maf1_feat = out_Features_cell(5).out_Features;

Msn2_5cond_time = [GD_Msn2_time, GD_Msn2_1_time, OS_Msn2_time, ZM_Msn2_time, PH_Msn2_time];
Msn2_5cond_track = [GD_Msn2_track, GD_Msn2_1_track, OS_Msn2_track, ZM_Msn2_track, PH_Msn2_track];
Msn4_5cond_time = [GD_Msn4_time, OS_Msn4_time, ZM_Msn4_time, PH_Msn4_time, ND_Msn4_time];
Msn4_5cond_track = [GD_Msn4_track, OS_Msn4_track, ZM_Msn4_track, PH_Msn4_track, ND_Msn4_track];
Dot6_5cond_time = [GD_Dot6_time, OS_Dot6_time, ZM_Dot6_time, PH_Dot6_time, ND_Dot6_time];
Dot6_5cond_track = [GD_Dot6_track, OS_Dot6_track, ZM_Dot6_track, PH_Dot6_track, ND_Dot6_track];
Maf1_5cond_time = [GD_Maf1_time, OS_Maf1_time, ZM_Maf1_time, PH_Maf1_time, ND_Maf1_time];
Maf1_5cond_track = [GD_Maf1_track, OS_Maf1_track, ZM_Maf1_track, PH_Maf1_track, ND_Maf1_track];
Pho4_5cond_time = [GD_Pho4_time, OS_Pho4_time, ZM_Pho4_time, PH_Pho4_time, ND_Pho4_time];
Pho4_5cond_track = [GD_Pho4_track, OS_Pho4_track, ZM_Pho4_track, PH_Pho4_track, ND_Pho4_track];

%% Mean Time Traces - filtered - Msn2
figure(1)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

condNames = {'GD', 'GD1', 'OS', 'ZM', 'PH'};
%%%% specify the subfolders present in the experiment %%%%
phases = {'Pre', 'Post'};
%%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
channel = '';
%%%% specify color of lines in plot %%%%
color_val = 'b'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};
%%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
val_to_plot = 'nf';

Nwells = length(Msn2_5cond_time);

%num_subplots = ceil(sqrt(Nwells));

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:Nwells % loop through each well
    
    subplot(2,3,theWells);
            
    all_tracks = Msn2_5cond_track{theWells};
    all_times = Msn2_5cond_time{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues_syc(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(condNames{theWells})
        %%%% specify the x and y range of the plot %%%%
        ylim([1 5])
        xlim([0 100])
        %axis([0 130 1 6])
    end
    
end
suptitle('Msn2')
print([base_dir, 'Msn2_5cond_MeanTrace.eps'], '-depsc')

% Msn4
figure(2)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};
%%%% specify the subfolders present in the experiment %%%%
phases = {'Pre', 'Post'};
%%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
channel = '';
%%%% specify color of lines in plot %%%%
color_val = 'b'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};
%%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
val_to_plot = 'nf';

Nwells = length(Msn4_5cond_time);

%num_subplots = ceil(sqrt(Nwells));

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:Nwells % loop through each well
    
    subplot(2,3,theWells);
            
    all_tracks = Msn4_5cond_track{theWells};
    all_times = Msn4_5cond_time{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues_syc(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(condNames{theWells})
        %%%% specify the x and y range of the plot %%%%
        ylim([1 5])
        xlim([0 100])
        %axis([0 130 1 6])
    end
    
end
suptitle('Msn4')
print([base_dir, 'Msn4_5cond_MeanTrace.eps'], '-depsc')

% Dot6 
figure(3)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};
%%%% specify the subfolders present in the experiment %%%%
phases = {'Pre', 'Post'};
%%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
channel = '';
%%%% specify color of lines in plot %%%%
color_val = 'b'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};
%%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
val_to_plot = 'nf';

Nwells = length(Dot6_5cond_time);

%num_subplots = ceil(sqrt(Nwells));

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:Nwells % loop through each well
    
    subplot(2,3,theWells);
            
    all_tracks = Dot6_5cond_track{theWells};
    all_times = Dot6_5cond_time{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues_syc(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(condNames{theWells})
        %%%% specify the x and y range of the plot %%%%
        ylim([1 5])
        xlim([0 100])
        %axis([0 130 1 6])
    end
    
end
suptitle('Dot6')
print([base_dir, 'Dot6_5cond_MeanTrace.eps'], '-depsc')

% Maf1
figure(4)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};
%%%% specify the subfolders present in the experiment %%%%
phases = {'Pre', 'Post'};
%%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
channel = '';
%%%% specify color of lines in plot %%%%
color_val = 'b'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};
%%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
val_to_plot = 'nf';

Nwells = length(Maf1_5cond_time);

%num_subplots = ceil(sqrt(Nwells));

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:Nwells % loop through each well
    
    subplot(2,3,theWells);
            
    all_tracks = Maf1_5cond_track{theWells};
    all_times = Maf1_5cond_time{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues_syc(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(condNames{theWells})
        %%%% specify the x and y range of the plot %%%%
        ylim([1 5])
        xlim([0 100])
        %axis([0 130 1 6])
    end
    
end
suptitle('Maf1')
print([base_dir, 'Maf1_5cond_MeanTrace.eps'], '-depsc')

% Pho4
figure(5)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};
%%%% specify the subfolders present in the experiment %%%%
phases = {'Pre', 'Post'};
%%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
channel = '';
%%%% specify color of lines in plot %%%%
color_val = 'b'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};
%%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
val_to_plot = 'nf';

Nwells = length(Pho4_5cond_time);

%num_subplots = ceil(sqrt(Nwells));

% "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
for theWells = 1:Nwells % loop through each well
    
    subplot(2,3,theWells);
            
    all_tracks = Pho4_5cond_track{theWells};
    all_times = Pho4_5cond_time{theWells};
    for pH = 1:length(phases) % loop through pre and post phases
        tracks = all_tracks.(phases{pH}); length(tracks)
        timevals = all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues_syc(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        title(condNames{theWells})
        %%%% specify the x and y range of the plot %%%%
        ylim([1 8])
        xlim([0 100])
        %axis([0 130 1 6])
    end
    
end
suptitle('Pho4')
print([base_dir, 'Pho4_5cond_MeanTrace.eps'], '-depsc')

%print('/Users/susanychen/Downloads/TempMatlabData20150808/MeanCellPlot_20150710PhosphateDeplet.eps', '-depsc')
%print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MeanCellPlot.pdf','-dpdf')

%% Concatenating structure
Msn2_5cond_feat(1).Msn2 = GD_Msn2_feat;
Msn2_5cond_feat(2).Msn2 = GD_Msn2_1_feat;
Msn2_5cond_feat(3).Msn2 = OS_Msn2_feat;
Msn2_5cond_feat(4).Msn2 = ZM_Msn2_feat;
Msn2_5cond_feat(5).Msn2 = PH_Msn2_feat;
Msn4_5cond_feat(1).Msn4 = GD_Msn4_feat;
Msn4_5cond_feat(2).Msn4 = OS_Msn4_feat;
Msn4_5cond_feat(3).Msn4 = ZM_Msn4_feat;
Msn4_5cond_feat(4).Msn4 = PH_Msn4_feat;
Msn4_5cond_feat(5).Msn4 = ND_Msn4_feat;
Dot6_5cond_feat(1).Dot6 = GD_Dot6_feat;
Dot6_5cond_feat(2).Dot6 = OS_Dot6_feat;
Dot6_5cond_feat(3).Dot6 = ZM_Dot6_feat;
Dot6_5cond_feat(4).Dot6 = PH_Dot6_feat;
Dot6_5cond_feat(5).Dot6 = ND_Dot6_feat;
Maf1_5cond_feat(1).Maf1 = GD_Maf1_feat;
Maf1_5cond_feat(2).Maf1 = OS_Maf1_feat;
Maf1_5cond_feat(3).Maf1 = ZM_Maf1_feat;
Maf1_5cond_feat(4).Maf1 = PH_Maf1_feat;
Maf1_5cond_feat(5).Maf1 = ND_Maf1_feat;
Pho4_5cond_feat(1).Pho4 = GD_Pho4_feat;
Pho4_5cond_feat(2).Pho4 = OS_Pho4_feat;
Pho4_5cond_feat(3).Pho4 = ZM_Pho4_feat;
Pho4_5cond_feat(4).Pho4 = PH_Pho4_feat;
Pho4_5cond_feat(5).Pho4 = ND_Pho4_feat;


%% Plotting "Number of Peaks" Feature

% Msn2
condNames1 = {'GD', 'GD1', 'OS', 'ZM', 'PH'};

numWells = length(Msn2_5cond_feat);
splot = ceil(sqrt(numWells));

clear uniqueVals
clear peakNumVecCell

for jj = 1:numWells
    
    display(strcat('This is jj: ', num2str(jj)))
    
    numCells = length(Msn2_5cond_feat(jj).Msn2);
    
    for ii = 1:numCells
        numPeaksVec(ii) = Msn2_5cond_feat(jj).Msn2(ii).cell.NumPeaks;
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
figure(6)
bar(1:numWells, initMatrix', 0.5,'stack') %[initMatrix(1,:)', initMatrix(2,:)', initMatrix(3,:)', initMatrix(4,:)'], 0.5, 'stack')
set(gca, 'xticklabel', condNames1)
xlabel('TRs'); ylabel('Number of Cells (count)')
legend(strsplit(num2str(diffPeaks)))
title('Proportion of Number of Peaks Msn2')
print([base_dir, 'Msn2_5cond_NumPeaks.eps'], '-depsc')

%print('/Users/susanychen/Downloads/TempMatlabData20150808/NumPeaksPlot_20150710PhosphateDeplet.eps', '-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Msn4
condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};

numWells = length(Msn4_5cond_feat);
splot = ceil(sqrt(numWells));

clear uniqueVals
clear peakNumVecCell

for jj = 1:numWells
    
    display(strcat('This is jj: ', num2str(jj)))
    
    numCells = length(Msn4_5cond_feat(jj).Msn4);
    
    for ii = 1:numCells
        numPeaksVec(ii) = Msn4_5cond_feat(jj).Msn4(ii).cell.NumPeaks;
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
figure(7)
bar(1:numWells, initMatrix', 0.5,'stack') %[initMatrix(1,:)', initMatrix(2,:)', initMatrix(3,:)', initMatrix(4,:)'], 0.5, 'stack')
set(gca, 'xticklabel', condNames)
xlabel('TRs'); ylabel('Number of Cells (count)')
legend(strsplit(num2str(diffPeaks)))
title('Proportion of Number of Peaks Msn4')
print([base_dir, 'Msn4_5cond_NumPeaks.eps'], '-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dot6
condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};

numWells = length(Dot6_5cond_feat);
splot = ceil(sqrt(numWells));

clear uniqueVals
clear peakNumVecCell

for jj = 1:numWells
    
    display(strcat('This is jj: ', num2str(jj)))
    
    numCells = length(Dot6_5cond_feat(jj).Dot6);
    
    for ii = 1:numCells
        numPeaksVec(ii) = Dot6_5cond_feat(jj).Dot6(ii).cell.NumPeaks;
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
figure(8)
bar(1:numWells, initMatrix', 0.5,'stack') %[initMatrix(1,:)', initMatrix(2,:)', initMatrix(3,:)', initMatrix(4,:)'], 0.5, 'stack')
set(gca, 'xticklabel', condNames)
xlabel('TRs'); ylabel('Number of Cells (count)')
legend(strsplit(num2str(diffPeaks)))
title('Proportion of Number of Peaks Dot6')
print([base_dir, 'Dot6_5cond_NumPeaks.eps'], '-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maf1
condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};

numWells = length(Maf1_5cond_feat);
splot = ceil(sqrt(numWells));

clear uniqueVals
clear peakNumVecCell

for jj = 1:numWells
    
    display(strcat('This is jj: ', num2str(jj)))
    
    numCells = length(Maf1_5cond_feat(jj).Maf1);
    
    for ii = 1:numCells
        numPeaksVec(ii) = Maf1_5cond_feat(jj).Maf1(ii).cell.NumPeaks;
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
figure(9)
bar(1:numWells, initMatrix', 0.5,'stack') %[initMatrix(1,:)', initMatrix(2,:)', initMatrix(3,:)', initMatrix(4,:)'], 0.5, 'stack')
set(gca, 'xticklabel', condNames)
xlabel('TRs'); ylabel('Number of Cells (count)')
legend(strsplit(num2str(diffPeaks)))
title('Proportion of Number of Peaks Maf1')
print([base_dir, 'Maf1_5cond_NumPeaks.eps'], '-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pho4
condNames = {'GD', 'OS', 'ZM', 'PH', 'ND'};

numWells = length(Pho4_5cond_feat);
splot = ceil(sqrt(numWells));

clear uniqueVals
clear peakNumVecCell

for jj = 1:numWells
    
    display(strcat('This is jj: ', num2str(jj)))
    
    numCells = length(Pho4_5cond_feat(jj).Pho4);
    
    for ii = 1:numCells
        numPeaksVec(ii) = Pho4_5cond_feat(jj).Pho4(ii).cell.NumPeaks;
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
figure(10)
bar(1:numWells, initMatrix', 0.5,'stack') %[initMatrix(1,:)', initMatrix(2,:)', initMatrix(3,:)', initMatrix(4,:)'], 0.5, 'stack')
set(gca, 'xticklabel', condNames)
xlabel('TRs'); ylabel('Number of Cells (count)')
legend(strsplit(num2str(diffPeaks)))
title('Proportion of Number of Peaks Pho4')
print([base_dir, 'Pho4_5cond_NumPeaks.eps'], '-depsc')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST PEAK ONLY:
% SUBPLOT: MAX HEIGHT, TIME MAX HEIGHT, PULSE WIDTH, OFFSLOPE, ONSLOPE

% Msn2
condNames = {'GD', 'GD1', 'OS', 'ZM' ,'PH'};

numWells = length(Msn2_5cond_feat);

first = 1;
clear maxHeightVecCell
clear maxHeightTimeVecCell
clear pulseWidthVecCell
clear offSlopeVecCell
clear onSlopeVecCell
for jj = 1:numWells
    clear maxHeightVec
    clear maxHeightTimeVec
    clear pulseWidthVec
    clear offSlopeVec
    clear onSlopeVec
    pp =1;
    
    numCells = length(Msn2_5cond_feat(jj).Msn2);
    for ii = 1:numCells
        %ii
        % filter out the no peak data
        if Msn2_5cond_feat(jj).Msn2(ii).cell.NumPeaks ~= 0
           
            % only find the first peak (could have more than 1 peak)!
            if Msn2_5cond_feat(jj).Msn2(ii).cell.NumPeaks >= 1 & ~isempty(Msn2_5cond_feat(jj).Msn2(ii).cell.MaxHeight(1))
                currMsn2 = Msn2_5cond_feat(jj).Msn2;
                maxHeightVec(pp) = currMsn2(ii).cell.MaxHeight(1);
                maxHeightTimeVec(pp) = currMsn2(ii).cell.TimeMaxHeight(1);
                pulseWidthVec(pp) = currMsn2(ii).cell.PulseWidth(1);
                offSlopeVec(pp) = currMsn2(ii).cell.OffSlope(1);
                onSlopeVec(pp) = currMsn2(ii).cell.OnSlope(1);
                
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
    
    maxHeightVecCell{jj} = maxHeightVec;
    maxHeightTimeVecCell{jj} = maxHeightTimeVec;
    pulseWidthVecCell{jj} = pulseWidthVec;
    offSlopeVecCell{jj} = offSlopeVec;
    onSlopeVecCell{jj} = onSlopeVec;
end

figure(11); 
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',condNames, 'yLabel', ' (nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('On Slope');
suptitle('Msn2 First Peak Features');

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 10 10])
print([base_dir, 'Msn2_5cond_firstPeaksFeats.eps'], '-depsc')

maxHeightVecCell_1 = maxHeightVecCell;
maxHeightTimeVecCell_1 = maxHeightTimeVecCell;
pulseWidthVecCell_1 = pulseWidthVecCell;
offSlopeVecCell_1 = offSlopeVecCell;
onSlopeVecCell_1 = onSlopeVecCell;


% Msn4
condNames = {'GD', 'OS', 'ZM', 'PH' ,'ND'};

numWells = length(Msn4_5cond_feat);

first = 1;
clear maxHeightVecCell
clear maxHeightTimeVecCell
clear pulseWidthVecCell
clear offSlopeVecCell
clear onSlopeVecCell

for jj = 1:numWells
    clear maxHeightVec
    clear maxHeightTimeVec
    clear pulseWidthVec
    clear offSlopeVec
    clear onSlopeVec
    pp =1;
    
    numCells = length(Msn4_5cond_feat(jj).Msn4);
    for ii = 1:numCells
        ii
        % filter out the no peak data
        if Msn4_5cond_feat(jj).Msn4(ii).cell.NumPeaks ~= 0
           
            % only find the first peak
            if Msn4_5cond_feat(jj).Msn4(ii).cell.NumPeaks >= 1 & first == 1
                currMsn4 = Msn4_5cond_feat(jj).Msn4;
                maxHeightVec(pp) = currMsn4(ii).cell.MaxHeight(1);
                maxHeightTimeVec(pp) = currMsn4(ii).cell.TimeMaxHeight(1);
                pulseWidthVec(pp) = currMsn4(ii).cell.PulseWidth(1);
                offSlopeVec(pp) = currMsn4(ii).cell.OffSlope(1);
                onSlopeVec(pp) = currMsn4(ii).cell.OnSlope(1);
                
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
    
    maxHeightVecCell{jj} = maxHeightVec;
    maxHeightTimeVecCell{jj} = maxHeightTimeVec;
    pulseWidthVecCell{jj} = pulseWidthVec;
    offSlopeVecCell{jj} = offSlopeVec;
    onSlopeVecCell{jj} = onSlopeVec;
end

figure(12); 
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', ' (nuc/cyt per 2min)'); title('On Slope');
suptitle('Msn4 First Peak Features');

maxHeightVecCell_1 = maxHeightVecCell;
maxHeightTimeVecCell_1 = maxHeightTimeVecCell;
pulseWidthVecCell_1 = pulseWidthVecCell;
offSlopeVecCell_1 = offSlopeVecCell;
onSlopeVecCell_1 = onSlopeVecCell;

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 10 10])
print([base_dir, 'Msn4_5cond_firstPeaksFeats.eps'], '-depsc')


% Dot6
condNames = {'GD', 'OS', 'ZM', 'PH' ,'ND'};

numWells = length(Dot6_5cond_feat);


first = 1;
clear maxHeightVecCell
clear maxHeightTimeVecCell
clear pulseWidthVecCell
clear offSlopeVecCell
clear onSlopeVecCell

for jj = 1:numWells
    clear maxHeightVec
    clear maxHeightTimeVec
    clear pulseWidthVec
    clear offSlopeVec
    clear onSlopeVec
    pp =1;
    
    numCells = length(Dot6_5cond_feat(jj).Dot6);
    for ii = 1:numCells
        ii
        % filter out the no peak data
        if Dot6_5cond_feat(jj).Dot6(ii).cell.NumPeaks ~= 0
           
            % only find the first peak
            if Dot6_5cond_feat(jj).Dot6(ii).cell.NumPeaks >= 1 & first == 1
                currDot6 = Dot6_5cond_feat(jj).Dot6;
                maxHeightVec(pp) = currDot6(ii).cell.MaxHeight(1);
                maxHeightTimeVec(pp) = currDot6(ii).cell.TimeMaxHeight(1);
                pulseWidthVec(pp) = currDot6(ii).cell.PulseWidth(1);
                offSlopeVec(pp) = currDot6(ii).cell.OffSlope(1);
                onSlopeVec(pp) = currDot6(ii).cell.OnSlope(1);
                
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
    
    maxHeightVecCell{jj} = maxHeightVec;
    maxHeightTimeVecCell{jj} = maxHeightTimeVec;
    pulseWidthVecCell{jj} = pulseWidthVec;
    offSlopeVecCell{jj} = offSlopeVec;
    onSlopeVecCell{jj} = onSlopeVec;
end

figure(13); 
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', ' (min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('On Slope');
suptitle('Dot6 First Peak Features');

maxHeightVecCell_1 = maxHeightVecCell;
maxHeightTimeVecCell_1 = maxHeightTimeVecCell;
pulseWidthVecCell_1 = pulseWidthVecCell;
offSlopeVecCell_1 = offSlopeVecCell;
onSlopeVecCell_1 = onSlopeVecCell;

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 10 10])
print([base_dir, 'Dot6_5cond_firstPeaksFeats.eps'], '-depsc')

% Maf1
condNames = {'GD', 'OS', 'ZM', 'PH' ,'ND'};

numWells = length(Maf1_5cond_feat);

first = 1;
clear maxHeightVecCell
clear maxHeightTimeVecCell
clear pulseWidthVecCell
clear offSlopeVecCell
clear onSlopeVecCell

for jj = 1:numWells
    clear maxHeightVec
    clear maxHeightTimeVec
    clear pulseWidthVec
    clear offSlopeVec
    clear onSlopeVec
    pp =1;
    
    numCells = length(Maf1_5cond_feat(jj).Maf1);
    for ii = 1:numCells
        ii
        % filter out the no peak data
        if Maf1_5cond_feat(jj).Maf1(ii).cell.NumPeaks ~= 0
           
            % only find the first peak
            if Maf1_5cond_feat(jj).Maf1(ii).cell.NumPeaks >= 1 & first == 1
                currMaf1 = Maf1_5cond_feat(jj).Maf1;
                maxHeightVec(pp) = currMaf1(ii).cell.MaxHeight(1);
                maxHeightTimeVec(pp) = currMaf1(ii).cell.TimeMaxHeight(1);
                pulseWidthVec(pp) = currMaf1(ii).cell.PulseWidth(1);
                offSlopeVec(pp) = currMaf1(ii).cell.OffSlope(1);
                onSlopeVec(pp) = currMaf1(ii).cell.OnSlope(1);
                
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
    
    maxHeightVecCell{jj} = maxHeightVec;
    maxHeightTimeVecCell{jj} = maxHeightTimeVec;
    pulseWidthVecCell{jj} = pulseWidthVec;
    offSlopeVecCell{jj} = offSlopeVec;
    onSlopeVecCell{jj} = onSlopeVec;
end

figure(14); 
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',condNames, 'yLabel', ' (nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('On Slope');
suptitle('Maf1 First Peak Features');

maxHeightVecCell_1 = maxHeightVecCell;
maxHeightTimeVecCell_1 = maxHeightTimeVecCell;
pulseWidthVecCell_1 = pulseWidthVecCell;
offSlopeVecCell_1 = offSlopeVecCell;
onSlopeVecCell_1 = onSlopeVecCell;

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 10 10])
print([base_dir, 'Maf1_5cond_firstPeaksFeats.eps'], '-depsc')

% Pho4
condNames = {'GD', 'OS', 'ZM', 'PH' ,'ND'};

numWells = length(Pho4_5cond_feat);

first = 1;
clear maxHeightVecCell
clear maxHeightTimeVecCell
clear pulseWidthVecCell
clear offSlopeVecCell
clear onSlopeVecCell

for jj = 1:numWells
    clear maxHeightVec
    clear maxHeightTimeVec
    clear pulseWidthVec
    clear offSlopeVec
    clear onSlopeVec
    pp =1;
    
    numCells = length(Pho4_5cond_feat(jj).Pho4);
    for ii = 1:numCells
        
        % filter out the no peak data
        if Pho4_5cond_feat(jj).Pho4(ii).cell.NumPeaks ~= 0
        
                figure;
                plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
                hold on;
                plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
                    
            % only find the first peak
            if Pho4_5cond_feat(jj).Pho4(ii).cell.NumPeaks >= 1 & first == 1
                currPho4 = Pho4_5cond_feat(jj).Pho4;
                maxHeightVec(pp) = currPho4(ii).cell.MaxHeight(1);
                maxHeightTimeVec(pp) = currPho4(ii).cell.TimeMaxHeight(1);
                pulseWidthVec(pp) = currPho4(ii).cell.PulseWidth(1);
                offSlopeVec(pp) = currPho4(ii).cell.OffSlope(1);
                onSlopeVec(pp) = currPho4(ii).cell.OnSlope(1);
                
                %if ii >30 & ii <32;
%                     figure;
%                     plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
%                     hold on;
%                     plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
                %end
                
                pp=pp+1;
            end
        end

    end
    
    maxHeightVecCell{jj} = maxHeightVec;
    maxHeightTimeVecCell{jj} = maxHeightTimeVec;
    pulseWidthVecCell{jj} = pulseWidthVec;
    offSlopeVecCell{jj} = offSlopeVec;
    onSlopeVecCell{jj} = onSlopeVec;
end

figure(15); 
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',condNames, 'yLabel', '(nuc/cyt per 2min)'); title('On Slope');
suptitle('Pho4 First Peak Features');

maxHeightVecCell_1 = maxHeightVecCell;
maxHeightTimeVecCell_1 = maxHeightTimeVecCell;
pulseWidthVecCell_1 = pulseWidthVecCell;
offSlopeVecCell_1 = offSlopeVecCell;
onSlopeVecCell_1 = onSlopeVecCell;

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 10 10])
print([base_dir, 'Pho4_5cond_firstPeaksFeats.eps'], '-depsc')