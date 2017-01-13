% set paths
clear
profile off
ipdir = '/Users/susanychen/GITREPOS/image_analysis2/';
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%% test how well the cell were processed/found
im =imread('RFP_p1_t1.tiff');
figure(1); imshow(im,[]); hold on;
for i=1:39;
plot(timecoursedata(1).celldata(i).Cyloc.RFP, timecoursedata(1).celldata(i).Cxloc.RFP, 'g.', 'markersize', 20); hold on;
end

%%
%figure(1)
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
fname_load = 'bPACpaperControl_09082012b.mat';

load([base_dir,fname_load]);%,'all_times_vec', 'all_tracks_vec')%,'all_times_vec','all_tracks_vec')
%%%% specify the subfolders present in the experiment %%%%
%%
% If "timecoursedata", then plot the average plot this way

%%%% specify color of lines in plot %%%%
color_val = 'g'; %cmap(jj,:);
%%%% specify plot line properties %%%%
plot_params = {'linewidth',1.5,'LineStyle','-'};

timevals = linspace(0,40, length(timecoursedata));

% plotting 
timetrace = [];
stdtrace = [];
for i = 1:length(timecoursedata)
meanvec = [];
for j=1:length(timecoursedata(i).celldata)
    %i
    %j
    meanvec(j) = timecoursedata(i).celldata(j).nf.RFP;
end
timetrace(i) = nanmedian(meanvec);
stdtrace(i) = nanstd(meanvec);
end

figure(2); errorbar(timevals,timetrace,stdtrace,'Color',color_val,plot_params{:}); 
hold on;
%plot(timevals, [4*ones(1,60), 5*ones(1,30), 4*ones(1,15), 5*ones(1,15)]);
plot(timevals, [4*ones(1,60), 5*ones(1,30), 4*ones(1,15), 4*ones(1,15)]);
axis([0 35 1 5])
%xlim = [0 150];
%ylim = [0 5];
phases = {'Estrodiol Induced Ras2 S24N @ T=0, Blue light betwn 20-30min'};
xlabel('Time (minutes)'); ylabel('Nuclear Localization (nuc/cyt)');
title(phases);

%print(strcat(base_dir, fname_load, '.pdf'),'-dpdf')

%%
%phases = {'Pre', 'Post'};
%phases = {'BL3min40uW'};
%phases = {'BL3min40uWcon'};
%phases = {'BL20min40uW'};
%phases = {'BL20min40uWcon'};
phases = {'ImageFiles'};

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
        timevals = linspace(0,40,120);%all_times.(phases{pH});

        [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
        hold on;
        %title(phases{pH})
        %%%% specify the x and y range of the plot %%%%
        %axis([0 130 1 6])
        ylim([1 5])
        xlim([0 40])
        xlabel('time (min)');
        ylabel('Msn2 nuclear localization (nuc/cyto)');
        title('Estrodiol induced RasDN @ T=0, Blue light btwn 20-30min, Blue light btwn 35-40min');
        plot(linspace(0,40,120), [4*ones(1,60), 5*ones(1,30), 4*ones(1,15), 5*ones(1,15)], 'b')
    end
    
end
print(strcat(base_dir, fname_load, '.pdf'),'-dpdf')
