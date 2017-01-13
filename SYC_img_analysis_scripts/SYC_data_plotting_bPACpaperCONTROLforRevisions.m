% set paths
clear
profile off
ipdir = '/Users/susanychen/GITREPOS/image_analysis2/';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%
figure(1)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = '/Users/susanychen/temp_uscopy_Datafolder/';
%%%% specify the .mat file to load %%%%
%fname_load = 'bPACpapercontrol_3min40uWJSO8-43.mat';
%fname_load = 'bPACpapercontrol_3min40uWJSO11-38.mat;
%fname_load = 'bPACpapercontrol_20min40uWJSO8-43.mat';
%fname_load = 'bPACpapercontrol_20min40uWJSO11-38.mat';
%fname_load = 'bPACpapercontrol_50min40uWJSO11-38.mat';
fname_load = 'bPACpapercontrol_50min40uWJSO8-43.mat';

load([base_dir,fname_load]);%,'all_times_vec', 'all_tracks_vec')%,'all_times_vec','all_tracks_vec')
%%%% specify the subfolders present in the experiment %%%%
%%
% If "timecoursedata", then plot the average plot this way

%%%% specify color of lines in plot %%%%
color_val = 'g'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};

timevals = [1:1:length(timecoursedata)];

% plotting 
timetrace = [];
stdtrace = [];
for i = 1:length(timecoursedata)
meanvec = [];
for j=1:length(timecoursedata(i).celldata)
    meanvec(j) = timecoursedata(i).celldata(j).nf.RFP;
end
timetrace(i) = nanmedian(meanvec);
stdtrace(i) = nanstd(meanvec);
end

figure(2); errorbar(timevals,timetrace,stdtrace,'Color',color_val,plot_params{:});
axis([0 70 1 5])
%xlim = [0 70];
%ylim = [0 5];
phases = {'BL50min40uW'};
title(phases);

print(strcat(base_dir, fname_load, '.pdf'),'-dpdf')

%%
%phases = {'Pre', 'Post'};
%phases = {'BL3min40uW'};
%phases = {'BL3min40uWcon'};
%phases = {'BL20min40uW'};
%phases = {'BL20min40uWcon'};
%phases = {'BL50min40uWcon'};

%%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
channel = 'RFP';
%%%% specify color of lines in plot %%%%
color_val = 'g'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};
%%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
val_to_plot = 'nf';

Nwells = length(all_tracks_vec);

num_subplots = ceil(sqrt(Nwells));

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
        title(phases{pH})
        %%%% specify the x and y range of the plot %%%%
        %axis([0 130 1 6])
        ylim([1 5])
        xlim([0 70])
    end
    
end
print(strcat(base_dir, fname_load, '.pdf'),'-dpdf')
