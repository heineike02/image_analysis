% SCRIPT to generate test plots - TEST PLOTS

% run the find cells script for only 1 well -- this is done with
% SCY_img_analysis_20150820_Dot6_7StressInputs_1SDC.m for example

% plot the tracks on top of image - make sure reasonable cells were found
im = imread('RFP_p4_t1.tiff');
figure; imshow(im,[]); hold on; count=0; for i = 1:length(all_tracks); if length(all_tracks_vec{1}.Post(i).Cyloc)>0.8*40; plot(all_tracks_vec{1}.Post(i).Cyloc,all_tracks_vec{1}.Post(i).Cxloc, 'g'); hold on; count=count+1; end; end;
title(['Dot6: P4: ',num2str(count), ' cells: 80% max trace length'])

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
% Pho4
ipdir = 'C:\Users\susanychen\GitHub\image_analysis\';
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150820_TR_7StressInputs\';
phases = {'Post'};
fname_load = 'Pho4_7StressInputs';
well_order = {'GD 0.05%','1M Sorb','ND', 'ND+AAD', '200ug/ml Zymo','Pi-deplete', 'pH8.7','SDC'};
perturbation = {'Pho4 with 7 different stresses'};
perMaxLength = 0.8;
SYC_filter_extraction_plotting_SC_template_1(ipdir, base_dir, phases, fname_load, well_order, perturbation, perMaxLength);

% Dot6

% Maf1
%%
