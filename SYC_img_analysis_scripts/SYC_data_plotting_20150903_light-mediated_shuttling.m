% SCRIPT to generate test plots - TEST PLOTS

% run the find cells script for only 1 well -- this is done with
% SCY_img_analysis_20150820_Dot6_7StressInputs_1SDC.m for example

%plot the tracks on top of image - make sure reasonable cells were found
im = imread('RFP_p1_t20.tiff');
figure; imshow(im,[]); hold on; count=0; 
for i = 1:length(all_tracks_vec{1}.Two); 
    %if length(all_tracks_vec{1}.Two(i).Cyloc)>0.8*40; 
        plot(all_tracks_vec{1}.Two(i).Cyloc(1),all_tracks_vec{1}.Two(i).Cxloc(1), 'g.', 'markersize',5); 
        hold on; 
        count=count+1; 
    %end; 
end;
title(['SYC84 pAdh1-mCherry-NESLOVNLS: P1: Inten4: ',num2str(count), ' cells: 80% max trace length'])

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
% 1. Get rid of doubly counted cells %% wasn't successful
fields = {'One','Two','Three'};
overlapDist = 6;
clear diff_dist_array;

for hh = 1:length(fieldnames(all_tracks_vec{1})); % pretend there is only 1 set
    all_tracks = all_tracks_vec{1}.(fields{hh}); 
    % go through each cell, and find the difference between its center
    % point
    display(length(all_tracks))
    
    for dd = 1:length(all_tracks)
        dist_at_dd = sqrt( all_tracks(dd).Cyloc(1)^2 + all_tracks(dd).Cxloc(1)^2 );
        for ee = 1:length(all_tracks)
        %yDelta = abs( all_tracks(dd).Cyloc(1) - all_tracks(dd+1).Cyloc(1) )
        %xDelta = abs( all_tracks(dd).Cxloc(1) - all_tracks(dd+1).Cxloc(1) )
        %deltaDist = sqrt( xDelta^2 + yDelta^2 )
        %diff_dist_array{dd,dd+1} = deltaDist;
        dist_at_ee = sqrt( all_tracks(ee).Cyloc(1)^2 + all_tracks(ee).Cxloc(1)^2 );
        diff_dist = abs( dist_at_dd - dist_at_ee );
        diff_dist_array{dd,ee} = diff_dist;
        % if the distance is between 0 and 9 (variable parameter), then
        % delete track
        %if deltaDist < overlapDist
        %    all_tracks_nodups{1}.(fields{hh})(dd) = [];
        %end
            if diff_dist ~= 0 & diff_dist < overlapDist;
                %all_tracks(ee).nf = zeros(1,length(all_tracks(ee).nf));
                all_tracks(ee).nf.RFP = 0;
            end  
        end   
    end
    
    all_tracks_nodubs{hh} = all_tracks;
    
    diff_dist_array_array{hh} = diff_dist_array;
    clear store_dubs_dd; clear store_dubs_ee; clear diff_dist_array;
end

tempArray = diff_dist_array_array{3};
tempArray = cell2mat(tempArray);


%%
% 2. filter the traces %%% did not filter as well
ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150902-20150903_light-mediated_shuttling\';
phases = {'One','Two','Three'};
fname_load = 'SYC84_lightmediatedshuttling_20150903_satOD';
well_order = {'Light'};
perturbation = {'Light perturbation to assess nuclear translocation'};
perMaxLength = 0.8;

% set paths
%clear
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
Xper = 0.8;
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
        all_tracks_filt.(phases{jj}) = filtered_tracks_struct;
        all_times_filt.(phases{jj}) = timevals;
        
        visual check
        if jj == 2
           figure; histogram(track_lengths(indc_of_tracks));
        end  
        
    end
    
    all_tracks_filt_vec{ii} = all_tracks_filt;
    all_times_filt_vec{ii} = all_times_filt;
end

save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',... % temporary adjustable
    'well_order', 'perturbation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% SATURATION OD
% 3. plot plots
save_img_file(1).filename = 'SingleCellTraces';
save_img_file(2).filename = 'NumberOfPeaks';
save_img_file(3).filename = 'FirstPeakFeatures';
fname_load = 'SYC84_lightmediatedshuttling_20150903_satOD';
val_to_plot = 'nf';

% initialize
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150902-20150903_light-mediated_shuttling\';;
experi_name = 'SYC84_lightmediatedshuttling_20150903_satOD';
load([base_dir,experi_name])

phases = {'One','Two','Three'};

% SC.MSN2-RFP plot individual cells (filtered, not smoothed)
% - variance of the time trace (captures presence of responders and non and
% different types of profiles)
% Msn2
figure(1)
clf

hold on
%%%% specify channel name %%%%
channel = 'RFP'; %'RFP'

sample_order = {'pSYC84 10minBL10minOFF Inten4','pSYC84 10minBL10minOFF Inten10','pSYC84 NoLight'};
 
Nfolders = length(fieldnames(all_tracks_vec{1}));
num_subplots = ceil(sqrt(Nfolders));

defineLight(1).light = [3*ones(1,20), 2.5*ones(1,20)];
defineLight(2).light = [3*ones(1,20), 2.5*ones(1,20)];
defineLight(3).light = [2.5*ones(1,20)];
%defineLight(4).light = [2.5*ones(1,20)];
%defineLight(5).light = [2.5*ones(1,20)];

for jj = 1:Nfolders
    jj
    all_tracks = all_tracks_vec{1}.(phases{jj});
    all_times = all_times_vec{1}.(phases{jj});

    %subplot(num_subplots,num_subplots,jj);
    subplot(1,3,jj);
    
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
    
    plot(all_times, defineLight(jj).light, 'c-.', 'linewidth',4);
    
    totCellNum = length(all_tracks);
    text(14,1.2,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    axis([0,20,1,4])
    title(sample_order{jj})
end
suptitle(save_img_file(1).filename);
print([base_dir, save_img_file(1).filename, ' 20150903 SatOD','.eps'], '-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOGARITHMIC OD
%%
% 3. plot plots
save_img_file(1).filename = 'SingleCellTraces';
save_img_file(2).filename = 'NumberOfPeaks';
save_img_file(3).filename = 'FirstPeakFeatures';
fname_load = 'SYC84_lightmediatedshuttling_20150903_logOD';
val_to_plot = 'nf';

% initialize
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150902-20150903_light-mediated_shuttling\';;
experi_name = 'SYC84_lightmediatedshuttling_20150903_logOD';
load([base_dir,experi_name])

phases = {'One','Two','Three'};

% SC.MSN2-RFP plot individual cells (filtered, not smoothed)
% - variance of the time trace (captures presence of responders and non and
% different types of profiles)
% Msn2
figure(1)
clf

hold on
%%%% specify channel name %%%%
channel = 'RFP'; %'RFP'

sample_order = {'pSYC84 10minBL10minOFF Inten4','pSYC84 10minBL10minOFF Inten10','pSYC84 NoLight'};
 
Nfolders = length(fieldnames(all_tracks_vec{1}));
num_subplots = ceil(sqrt(Nfolders));

defineLight(1).light = [3*ones(1,20), 2.5*ones(1,20)];
defineLight(2).light = [3*ones(1,20), 2.5*ones(1,20)];
defineLight(3).light = [2.5*ones(1,20)];
%defineLight(4).light = [2.5*ones(1,20)];
%defineLight(5).light = [2.5*ones(1,20)];

for jj = 1:Nfolders
    jj
    all_tracks = all_tracks_vec{1}.(phases{jj});
    all_times = all_times_vec{1}.(phases{jj});

    %subplot(num_subplots,num_subplots,jj);
    subplot(1,3,jj);
    
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
    
    plot(all_times, defineLight(jj).light, 'c-.', 'linewidth',4);
    
    totCellNum = length(all_tracks);
    text(14,1.2,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    axis([0,20,1,4])
    title(sample_order{jj})
end
suptitle(save_img_file(1).filename);
print([base_dir, save_img_file(1).filename, ' 20150903 LogOD','.eps'], '-depsc')