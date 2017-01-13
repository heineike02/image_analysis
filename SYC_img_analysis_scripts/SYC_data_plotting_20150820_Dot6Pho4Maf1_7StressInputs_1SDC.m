% SCRIPT to generate test plots - TEST PLOTS

% run the find cells script for only 1 well -- this is done with
% SCY_img_analysis_20150820_Dot6_7StressInputs_1SDC.m for example

% plot the tracks on top of image - make sure reasonable cells were found
im = imread('RFP_p4_t1.tiff');
figure; imshow(im,[]); hold on; count=0; for i = 1:length(all_tracks); if length(all_tracks_vec{1}.Post(i).Cyloc)>0.8*40; plot(all_tracks_vec{1}.Post(i).Cyloc,all_tracks_vec{1}.Post(i).Cxloc, 'g'); hold on; count=count+1; end; end;
title(['Dot6: P4: ',num2str(count), ' cells: 80% max trace length'])

% plot smoothing using moving filter
% before executing this code, run the
% "SYC_filter_extraction_plotting_SC_template_1" code first
figure(1); filters = {'moving', 'lowess', 'sgolay'}; 
for i= 1:25
subplot(5,5,i); plot(all_tracks_filt_vec{1}.Post(i).nf.RFP,'b'); hold on; 
plot(smooth(all_tracks_filt_vec{1}.Post(i).nf.RFP, filters{1}), 'r'); hold on;
end
suptitle('Dot6: P1: Smoothing with Moving Average filter Span5')

% plot peak finding traces
figure(2); for i = 1:81; tempVar = out_Features_cell(1).out_Features(i).cell;
subplot(9,9,i); plot(tempVar.TimeTraceTimes, tempVar.TimeTrace); hold on; plot(tempVar.TimeMaxHeight, tempVar.MaxHeight, 'r.'); axis([0 120 1 6]);
end
suptitle('Dot6: P1: Peak Finding with Med+sel(0.1) threshold');

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

%% Notes: Before commiting this next round of scripts: remember to pull out any arbitrary parameters and explicitly define them in the script that runs the function
%% Notes: Important to keep improving the peak finding algorithm

%%
% 1. run raw images through the "image analysis template". Done in a
% different script
%%
% 2. filter the traces
% Dot6
ipdir = '/Users/susanychen/GITREPOS/image_analysis2/';
base_dir = '/Users/susanychen/20150820_TR_7StressInputs/';
phases = {'Post'};
fname_load = 'Dot6_7StressInputs';
well_order = {'GD 0.05%','1M Sorb','ND', 'ND+AAD', '200ug/ml Zymo','Pi-deplete', 'pH8.7','SDC'};
perturbation = {'Dot6 with 7 different stresses'};
perMaxLength = 0.8;
smoothFn.method = 'moving';
smoothFn.span = 30;
smoothFn.flag = 1;
SYC_filter_extraction_plotting_SC_template_1(ipdir, base_dir, phases, fname_load, well_order, perturbation, perMaxLength, smoothFn);

% Pho4
ipdir = '/Users/susanychen/Documents/image_analysis';
base_dir = '/Users/susanychen/Downloads/tempData20150820/';
phases = {'Post'};
fname_load = 'Pho4_7StressInputs';
well_order = {'GD 0.05%','1M Sorb','ND', 'ND+AAD', '200ug/ml Zymo','Pi-deplete', 'pH8.7','SDC'};
perturbation = {'Pho4 with 7 different stresses'};
perMaxLength = 0.8;
smoothFn.method = 'moving';
smoothFn.span = 30;
smoothFn.flag = 1;
SYC_filter_extraction_plotting_SC_template_1(ipdir, base_dir, phases, fname_load, well_order, perturbation, perMaxLength, smoothFn);

% Maf1
ipdir = '/Users/susanychen/Documents/image_analysis';
base_dir = '/Users/susanychen/Downloads/tempData20150820/';
phases = {'Post'};
fname_load = 'Maf1_7StressInputs';
well_order = {'GD 0.05%','1M Sorb','ND', 'ND+AAD', '200ug/ml Zymo','Pi-deplete', 'pH8.7','SDC'};
perturbation = {'Maf1 with 7 different stresses'};
perMaxLength = 0.8;
smoothFn.method = 'moving';
smoothFn.span = 30;
smoothFn.flag = 1;
SYC_filter_extraction_plotting_SC_template_1(ipdir, base_dir, phases, fname_load, well_order, perturbation, perMaxLength, smoothFn);
%%
%% Plot raw mean traces
%% SC.Msn2-RFP Mean Time Traces

addpath('/Users/susanychen/GITREPOS/image_analysis2')

% Dot6
figure(1)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'GD 0.05%','1M Sorb','ND','ND+AAD','200ug/ml Zymo',...
    'Pi-deplete','pH8.7'};

base_dir = '/Users/susanychen/20150820_TR_7StressInputs/';
fname_save = 'Dot6_7StressInputs.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [ 0,0,1;  %Blue
    1,0,0; %Red
    0,1,0;
    1,1,0;
    0,0,0;
    0,1,1;
    1,0.5,1];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
phases = {'Post'};
perm = [1,2,3,4,5,6,7];
N_RFP = length(perm)
legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_RFP_plot),1);
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
%htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('Dot6 Localization Stress Panel')
xlabel('Time (min)')
ylabel('Nuclear Localization (nuc/cyt)')
axis([0 120 1 8])

% Maf1
figure(2)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'GD 0.05%','1M Sorb','ND','ND+AAD','200ug/ml Zymo',...
    'Pi-deplete','pH8.7'};

base_dir = '/Users/susanychen/20150820_TR_7StressInputs/';
fname_save = 'Maf1_7StressInputs.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [ 0,0,1;  %Blue
    1,0,0; %Red
    0,1,0;
    1,1,0;
    0,0,0;
    0,1,1;
    1,0.5,1];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
phases = {'Post'};
perm = [1,2,3,4,5,6,7];
N_RFP = length(perm)
legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_RFP_plot),1);
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
%htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('Maf1 Localization Stress Panel')
xlabel('Time (min)')
ylabel('Nuclear Localization (nuc/cyt)')
axis([0 120 1 8])

% Pho4
figure(3)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'GD 0.05%','1M Sorb','ND','ND+AAD','200ug/ml Zymo',...
    'Pi-deplete','pH8.7'};

base_dir = '/Users/susanychen/20150820_TR_7StressInputs/';
fname_save = 'Pho4_7StressInputs.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [ 0,0,1;  %Blue
    1,0,0; %Red
    0,1,0;
    1,1,0;
    0,0,0;
    0,1,1;
    1,0.5,1];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
phases = {'Post'};
perm = [1,2,3,4,5,6,7];
N_RFP = length(perm)
legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_RFP_plot),1);
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
%htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('Pho4 Localization Stress Panel')
xlabel('Time (min)')
ylabel('Nuclear Localization (nuc/cyt)')
axis([0 120 1 8])

%%
% 3. 4. 5. 6. plot the different traces/characteristics
save_img_file(1).filename = 'Dot6 Single Cell Traces';
save_img_file(2).filename = 'Dot6 Proportion Peak Number';
save_img_file(3).filename = 'Dot6 First Peak Features';
fname_load = 'Dot6_7StressInputs_filter';
plotVars.axis = [0 120 0 6.5];
SYC_filter_extraction_plotting_SC_template_2(base_dir, fname_load, well_order, phases, save_img_file,plotVars)

% Pho4
save_img_file(1).filename = 'Pho4 Single Cell Traces';
save_img_file(2).filename = 'Pho4 Proportion Peak Number';
save_img_file(3).filename = 'Pho4 First Peak Features';
fname_load = 'Pho4_7StressInputs_filter';
plotVars.axis = [0 120 0 15];
SYC_filter_extraction_plotting_SC_template_2(base_dir, fname_load, well_order, phases, save_img_file,plotVars)

% Maf1
save_img_file(1).filename = 'Maf1 Single Cell Traces';
save_img_file(2).filename = 'Maf1 Proportion Peak Number';
save_img_file(3).filename = 'Maf1 First Peak Features';
fname_load = 'Maf1_7StressInputs_filter';
plotVars.axis = [0 120 0 12];
SYC_filter_extraction_plotting_SC_template_2(base_dir, fname_load, well_order, phases, save_img_file, plotVars)