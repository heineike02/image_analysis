% SCRIPT to generate test plots - TEST PLOTS

% run the find cells script for only 1 well -- this is done with
% SCY_img_analysis_20150820_Dot6_7StressInputs_1SDC.m for example

% %plot the tracks on top of image - make sure reasonable cells were found
im = imread('RFP_p1_t24.tiff');
figure; imshow(im,[]);
hold on;
for i = 1:length(all_tracks_vec{1}.Dot6)
plot(all_tracks_vec{1}.Dot6(i).Cyloc(1),all_tracks_vec{1}.Dot6(i).Cxloc(1), 'g.'); hold on;
end

% % the other plot track
% first load the .mat file
im = imread('RFP_p1_t23.tiff');
figure; imshow(im,[]);
hold on;
for i = 1:length(all_tracks_vec{1}.Tod6Cbk1)
plot(all_tracks_vec{1}.Tod6Cbk1(i).Cyloc(1),all_tracks_vec{1}.Tod6Cbk1(i).Cxloc(1), 'g.'); hold on;
end
% plot peak finding traces

%%
% README - instructions on how to process image data
% 1. run raw images through the "image analysis template" (make sure
% reasonable cells were found) - finds the cut-off threshold for traces
% that are real cells (same filtering for length)
% 2. filter the traces to only cells with longer traces
% 3. plot the median with std
% 4. plot the single cell traces (include trace numbers)
% 5. plot the proportion of peaks found
% 6. plot the 1st peak characteristics
%%
% 1. run raw images through the "image analysis template". Done in a
% different script
%%
% 2. filter the traces

%%
% 3. plot plots
% 3. plot plots
save_img_file(1).filename = 'SingleCellTraces';
save_img_file(2).filename = 'NumberOfPeaks';
save_img_file(3).filename = 'FirstPeakFeatures';
fname_load = 'LOVnoLOVlight_20150929';
val_to_plot = 'nf';
%val_to_plot = 'nmi';

% initialize
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150929_TR-NLSvar11_LOVnoLOV_light_GD_experiment\';
experi_name = fname_load;
load([base_dir,experi_name])

phases =  {'one','two','three','four', 'five', 'eight', 'nine'};

well_order = {'Dot6-mCherry-NLSvar11 BMH35', 'Tod6-mCherry-NLSvar11 BMH35', 'Pho4-mCherry-NLSvar11 BMH35', 'Dot6-mCherry BMH35', 'Pho4-mCherry BMH35',  'Msn4-mCherry BMH35','Dot6-mCherry-NLSvar11 BMH35 NO LIGHT'};
sample_order = well_order;

% SC.MSN2-RFP plot individual cells (filtered, not smoothed)
% - variance of the time trace (captures presence of responders and non and
% different types of profiles)
% Msn2
figure(1)
clf

hold on
%%%% specify channel name %%%%
channel = 'RFP'; %'RFP'
 
Nfolders = length(fieldnames(all_tracks_vec{1}));
num_subplots = ceil(sqrt(Nfolders));

defineLight(1).light = [2*ones(1,4),3*ones(1,10), 2*ones(1,10)];

%defineLight(1).light = [200*ones(1,4),500*ones(1,10), 200*ones(1,10)];


%defineLight(1).light = [3*ones(1,6),4*ones(1,6), 3*ones(1,6), 4*ones(1,6), 3*ones(1,6)];

theAxis{1} = [0 12 1.5 18];
%theAxis{7} = [0 20 1.5 3];
%theAxis{8} = [0 12 1.5 3];
%theAxis{9} = [0 8 1.5 3];
%theAxis{10} = [0 4 1.5 3];
%theAxis{11} = [0 2 1.5 3];
for jj = 1:Nfolders
    jj
    all_tracks = all_tracks_vec{1}.(phases{jj});
    all_times = all_times_vec{1}.(phases{jj});

    %subplot(num_subplots,num_subplots,jj);
    subplot(3,3,jj);
    
    % loop through single cells
    hold all
    xlabel('time(min)')
    for nn = 1:length(all_tracks)
       if isempty(channel)
          y_vec = [all_tracks(nn).(val_to_plot)];
       else
          y_vec = [smooth(all_tracks(nn).(val_to_plot).(channel))];
       end
       time_inds = all_tracks(nn).times;
       t_vec = all_times(time_inds);
       fig_out = plot(t_vec,y_vec);
       %patchline
       hold on;
    end
    
    plot(all_times, defineLight(1).light, 'c-.', 'linewidth',4);
    
    totCellNum = length(all_tracks);
    text(1,1,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    %if jj <=6
        %axis(theAxis{1});
        xlim([0 12])
    %else
        %axis(theAxis{jj});
    %end
    
    title(sample_order{jj})
end
suptitle(save_img_file(1).filename);
%print([base_dir, save_img_file(1).filename, ' 20150929 TR-NLSvar11 LOVnoLOV light NF','.eps'], '-depsc')
%print([base_dir, save_img_file(1).filename, ' 20150929 TR-NLSvar11 LOVnoLOV light NMI','.eps'], '-depsc')

%% Plot the mean with standard deviation

figure(1)
clf

hold on
%%%% specify channel name %%%%
channel = 'RFP'; %'RFP'
cmap_RFP = [ 0,0,1;  %Blue
    1,0,0; %Red
    0 1 0; %Green
    0,0,1;  %Blue
    0,1,0; %Red
    0 1 0; %Green
    0,0,1;  %Blue
    ];  
plot_params = {'linewidth',1.5,'LineStyle','-'};

Nfolders = length(fieldnames(all_tracks_vec{1}));
num_subplots = ceil(sqrt(Nfolders));

defineLight(1).light = [2*ones(1,4),3*ones(1,10), 2*ones(1,10)];
%defineLight(1).light = [200*ones(1,4),500*ones(1,10), 200*ones(1,10)];
%defineLight(1).light = [3*ones(1,6),4*ones(1,6), 3*ones(1,6), 4*ones(1,6), 3*ones(1,6)];

theAxis{1} = [0 12 1.5 18];
%theAxis{7} = [0 20 1.5 3];
%theAxis{8} = [0 12 1.5 3];
%theAxis{9} = [0 8 1.5 3];
%theAxis{10} = [0 4 1.5 3];
%theAxis{11} = [0 2 1.5 3];
for jj = 1:Nfolders
    jj
    all_tracks = all_tracks_vec{1}.(phases{jj});
    all_times = all_times_vec{1}.(phases{jj});

    %subplot(num_subplots,num_subplots,jj);
    subplot(3,3,jj);
    
    % loop through single cells
    hold all
    xlabel('time(min)')
    xlim([0 12])
    ylabel('nuc enrich (nuc/cyt)')
%     for nn = 1:length(all_tracks)
%        if isempty(channel)
%           y_vec = [all_tracks(nn).(val_to_plot)];
%        else
%           y_vec = [smooth(all_tracks(nn).(val_to_plot).(channel))];
%        end
%        time_inds = all_tracks(nn).times;
%        t_vec = all_times(time_inds);
%        fig_out = plot(t_vec,y_vec);
%        %patchline
%        hold on;
%     end
    color_val = cmap_RFP(jj,:);
    
    p = plot_meanvalues_syc(all_times,all_tracks,channel,color_val,0,'nf','plot_params',plot_params)
    
    %plot(all_times, defineLight(1).light, 'c-.', 'linewidth',4);
    
    totCellNum = length(all_tracks);
    text(1,10,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    %if jj <=6
        %axis(theAxis{1});
        %xlim([0 12])
    %else
        %axis(theAxis{jj});
    %end
    
    title(sample_order{jj})
end
suptitle(save_img_file(1).filename);

%%
% 3. plot plots
% This is for the 2 other samples that are RFP and YFP
save_img_file(1).filename = 'SingleCellTraces';
save_img_file(2).filename = 'NumberOfPeaks';
save_img_file(3).filename = 'FirstPeakFeatures';
fname_load = 'LOVnoLOVGD_TR-NLSvar11_20150929';
val_to_plot = 'nf';
%val_to_plot = 'nmi';

% initialize
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150929_TR-NLSvar11_LOVnoLOV_light_GD_experiment\';
experi_name = fname_load;
load([base_dir,experi_name])

phases =  {'beforeGD','afterGD'};

well_order = {'Dot6-mCherry-NLSvar11 BMH35', 'Tod6-mCherry-NLSvar11 BMH35', 'Pho4-mCherry-NLSvar11 BMH35', 'Dot6-mCherry BMH35', 'Pho4-mCherry BMH35'};
sample_order = well_order;

% SC.MSN2-RFP plot individual cells (filtered, not smoothed)
% - variance of the time trace (captures presence of responders and non and
% different types of profiles)
% Msn2
figure(1)
clf

hold on
%%%% specify channel name %%%%
%channel = 'RFP'; %'RFP'
channel = 'RFP';
 
Nfolders = length(all_tracks_vec);
num_subplots = ceil(sqrt(Nfolders));

%defineLight(1).light = [2*ones(1,80)];%
%defineLight(1).light = [200*ones(1,4),300*ones(1,10), 200*ones(1,120)];


%defineLight(1).light = [3*ones(1,6),4*ones(1,6), 3*ones(1,6), 4*ones(1,6), 3*ones(1,6)];

theAxis{1} = [0 80 1.5 12];
%theAxis{7} = [0 20 1.5 3];
%theAxis{8} = [0 12 1.5 3];
%theAxis{9} = [0 8 1.5 3];
%theAxis{10} = [0 4 1.5 3];
%theAxis{11} = [0 2 1.5 3];
for jj = 1:Nfolders
    jj
    all_tracks = all_tracks_vec{jj}.(phases{2});
    all_times = all_times_vec{jj}.(phases{2});

    %subplot(num_subplots,num_subplots,jj);
    subplot(2,3,jj);
    
    % loop through single cells
    hold all
    xlabel('time(min)')
    for nn = 1:length(all_tracks)
       if isempty(channel)
          y_vec = [all_tracks(nn).(val_to_plot)];
       else
          y_vec = [smooth(all_tracks(nn).(val_to_plot).(channel))];
       end
       time_inds = all_tracks(nn).times;
       t_vec = all_times(time_inds);
       fig_out = plot(t_vec,y_vec);
       %patchline
       hold on;
    end
    
    %plot(all_times, defineLight(jj).light, 'c-.', 'linewidth',4);
    
    totCellNum = length(all_tracks);
    text(1,900,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    %if jj <=6
        %axis(theAxis{1});
    %else
    %    axis(theAxis{jj});
    %end
    
    title(sample_order{jj})
end
suptitle(save_img_file(1).filename);
%print([base_dir, save_img_file(1).filename, ' 20150929 TR-NLSvar11 LOVnoLOV GD NMI','.eps'], '-depsc')

%% plot the mean with standard deviation

figure(1)
clf

hold on
%%%% specify channel name %%%%
channel = 'RFP'; %'RFP'
cmap_RFP = [ 0,0,1;  %Blue
    1,0,0; %Red
    0 1 0; %Green
    0,0,1;  %Blue
    0,1,0; %Red
    0 1 0; %Green
    0,0,1;  %Blue
    ];  
plot_params = {'linewidth',1.5,'LineStyle','-'};

Nfolders = length(all_tracks_vec);
num_subplots = ceil(sqrt(Nfolders));

defineLight(1).light = [2*ones(1,4),3*ones(1,10), 2*ones(1,10)];
%defineLight(1).light = [200*ones(1,4),500*ones(1,10), 200*ones(1,10)];
%defineLight(1).light = [3*ones(1,6),4*ones(1,6), 3*ones(1,6), 4*ones(1,6), 3*ones(1,6)];

theAxis{1} = [0 12 1.5 18];
%theAxis{7} = [0 20 1.5 3];
%theAxis{8} = [0 12 1.5 3];
%theAxis{9} = [0 8 1.5 3];
%theAxis{10} = [0 4 1.5 3];
%theAxis{11} = [0 2 1.5 3];
for jj = 1:Nfolders
    jj
    all_tracks = all_tracks_vec{jj}.(phases{2});
    all_times = all_times_vec{jj}.(phases{2});

    %subplot(num_subplots,num_subplots,jj);
    subplot(3,3,jj);
    
    % loop through single cells
    hold all
    xlabel('time(min)')
   %xlim([0 12])
    ylabel('nuc enrich (nuc/cyt)')
%     for nn = 1:length(all_tracks)
%        if isempty(channel)
%           y_vec = [all_tracks(nn).(val_to_plot)];
%        else
%           y_vec = [smooth(all_tracks(nn).(val_to_plot).(channel))];
%        end
%        time_inds = all_tracks(nn).times;
%        t_vec = all_times(time_inds);
%        fig_out = plot(t_vec,y_vec);
%        %patchline
%        hold on;
%     end
    color_val = cmap_RFP(jj,:);
    
    p = plot_meanvalues_syc(all_times,all_tracks,channel,color_val,0,'nf','plot_params',plot_params)
    
    %plot(all_times, defineLight(1).light, 'c-.', 'linewidth',4);
    
    totCellNum = length(all_tracks);
    text(1,10,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    %if jj <=6
        %axis(theAxis{1});
        %xlim([0 12])
    %else
        %axis(theAxis{jj});
    %end
    
    title(sample_order{jj})
end
suptitle(save_img_file(1).filename);