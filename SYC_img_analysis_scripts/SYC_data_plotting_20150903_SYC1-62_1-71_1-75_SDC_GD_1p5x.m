% SCRIPT to generate test plots - TEST PLOTS

% run the find cells script for only 1 well -- this is done with
% SCY_img_analysis_20150820_Dot6_7StressInputs_1SDC.m for example

% %plot the tracks on top of image - make sure reasonable cells were found
% im = imread('RFP_p4_t1.tiff');
% figure; imshow(im,[]);
% hold on;
% for i = 1:length(all_tracks_vec{1}.One)
% plot(all_tracks_vec{1}.One(i).Cyloc(1),all_tracks_vec{1}.One(i).Cxloc(1), 'g.'); hold on;
% end
% title(['Dot6: P4: ',num2str(count), ' cells: 80% max trace length'])

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
ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150903_SYC1-62_1-71_1-75_SDC_GD_1p5x\';
phases = {'Post'};
fname_load = 'SYC1-62_1-71_1-75_SDC_GD_1p5x';
well_order = {'Dot6 SDC','Pho4 SDC','Maf1 SDC','Dot6 0.05%glucose','Pho4 0.05%glucose','Maf1 0.05%glucose'};
perturbation = {'Dot6, Pho4, Maf1 with SDC or 0.05% glucose'};
perMaxLength = 0.8;
SYC_filter_extraction_plotting_SC_template_1(ipdir, base_dir, phases, fname_load, well_order, perturbation, perMaxLength);

%%
% 3. plot plots
save_img_file(1).filename = 'SingleCellTraces';
save_img_file(2).filename = 'NumberOfPeaks';
save_img_file(3).filename = 'FirstPeakFeatures';
fname_load = 'SYC1-62_1-71_1-75_SDC_GD_1p5x_filter';

SYC_filter_extraction_plotting_SC_template_2(base_dir, fname_load, well_order, phases, save_img_file)