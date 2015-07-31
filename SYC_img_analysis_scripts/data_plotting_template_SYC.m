function all_tracks = data_plotting_template()
% Note that this function EXTRACTS SINGLE CELL FEATURES and PLOTS THE TIME
% TRACES AND THE FEATURES

% basic initial values
clear
profile off
ipdir = 'C:\Users\susanychen\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\'
%%%% phases = subfolders in the same experiment %%%%
phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
%%%% fname_load = the .mat folder that contains the single cell data %%%%
fname_load = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720.mat';

%% Filter single cells traces
% at least 50% of max length

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
        
%%%% fname_save = save these filtered traces %%%%
fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save],'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec')

%% Feature Extraction and (Optional) Smoothing

%%%% specify the right folder where the filtered traces live %%%%
fname_load = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
load([fname_load],'all_tracks_filt_vec', 'all_times_filt_vec','posvec')
% clear old variables
clear t_singlecell_cell; clear y_singlecell_cell;

% set variables
%%%% If only 1 channel, leave blank, if more than 1 channel, specify the
%%%% channel name %%%%
channel = ''; %'RFP'
Nwells = length(all_tracks_filt_vec);

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

% set the variables for feature extraction -- stopped here
var_to_plot = 'nf';
smoothFn.flag = 1;
smoothFn.method = 'sgolay';
smoothFn.span = 0.4; % span of 40%

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
        
        peakFindr.Sel = 0.08; % <-- this is the golden parameter for best peak detection %(max(singleCells1.nf) - min(singleCells1.nf))/8; %nanmean(singleCells1.nf)% + 0.125*nanstd(singleCells1.nf); % wiggle these 2 parameters
        peakFindr.Thresh = 1.5;
        peakFindr.StdMinus = 0; %0.2; % wiggle these 2 parameters
        
        %% test plot
        %figure; plot(singleCells1.nf);
        %%plot(out_Features(kk).cell.TimeMaxHeight', out_Features(kk).cell.MaxHeight, 'ro');
        
        % extract features: number of peaks, max height(s), time(s) of max height, rise time (on slope, off slope, on slope time, off slope time), duration of pulse
        out_Features(kk).cell = extract_singlecell_features_Analysis(tsingleCells, singleCells1,var_to_plot, smoothFn, peakFindr);
        
        %% test plot
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
    out_Features_cell(ii).out_Features = out_Features;
%     out_Features_cell(ii).var_over_time =1; 
%     out_Features_cell(ii).numPeaksVar =1;
%     out_Features_cell(ii).maxHeightVar = 1;
%     out_Features_cell(ii).timeMaxHeightVar = 1;
%     out_Features_cell(ii).onSlope
    clear out_Features;
end

%%%% TR_names =  the labels that go with each well %%%%
TR_names = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
    'SYC1-70 crz1', 'SYC1-62 Dot6', ...
    'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
    'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
    'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
    'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};
%%%% stress_type = labels of stress type %%%%
stress_type = {'20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic'...
    '20150710pHbasic','20150710pHbasic','20150710pHbasic'};

%%%% fname_save = save plotting values/features %%%%
fname_save = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
save([fname_save], 'all_tracks_vec', 'all_times_vec','all_tracks_filt_vec', 'all_times_filt_vec','posvec',...
    'out_Features_cell',...
    'TR_names','stress_type')

%% Mean Time Traces - filtered SYC (filtered, not smoothed)
figure(1)
clf
%hold on
clear mean_val_cell; clear std_val_cell;

%%%% base_dir = directory where single cell data is stored %%%%
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
%%%% specify the .mat file to load %%%%
fname_load = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720_filter.mat';
load([base_dir,fname_load],'all_times_filt_vec','all_tracks_filt_vec','out_Features_cell','TR_names', 'stress_type')
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

% % should be pre-specified
% TR_names = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
%     'SYC1-70 crz1', 'SYC1-62 Dot6', ...
%     'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
%     'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
%     'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
%     'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};

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
        %%%% specify the x and y range of the plot %%%%
        axis([0 130 1 6])
    end
    
end
print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MeanCellPlot.eps','-depsc')

%% SC.MSN2-RFP plot individual cells (filtered, not smoothed)
%% Note the variance of single cell traces
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

% % should be prespecified
% TR_names = {'SYC1-67 Cad1 phosphate deplete SDC', 'SYC1-73 Com2', ...
%     'SYC1-70 crz1', 'SYC1-62 Dot6', ...
%     'SYC1-75 Maf1', 'SYC1-72 Msn2', ...
%     'SYC1-64 Yap1', 'SYC1-76 Stb3', ...
%     'SYC1-69 Sko1', 'SYC1-65 Rtg3', ...
%     'SYC1-71 Pho4', 'SYC1-68 Nrg2','SYC1-74 Msn4'};

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
    title(TR_names{jj})
end
print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\SingleCellPlot.eps','-depsc')

% % % %% Cluster the single cell traces into a few groups
% % % Nwells_singleCell = y_singlecell_cell(:,2);
% % % Nwells_t_singleCell = t_singlecell_cell(:,2);
% % %     singleCell = Nwells_singleCell{2};
% % %     tsingleCell = Nwells_t_singleCell{2};
% % %     
% % %     
% % %     for tt = 1:length(singleCell)
% % %         [Dist(tt),D,k(tt),w]=dtw(singleCell(1).nf',singleCell(tt).nf');
% % %     end
% % %     
% % %     [idx,C]=kmeans(Dist./k,2);
% % %     ind1=find(idx==1);
% % %     ind2=find(idx==2);
% % %     ind3=find(idx==3);
% % %     ind4=find(idx==4);
% % %     
% % %     figure(1); hold on;
% % %     for yy = 1:length(ind1)
% % %         plot(singleCell(ind1(yy)).nf)
% % %     end
% % %     axis([0 60 1 7])
% % %     
% % %     figure(2); hold on;
% % %     for yy = 1:length(ind2)
% % %         plot(singleCell(ind2(yy)).nf)
% % %     end
% % %     axis([0 60 1 7])
% % % 
% % %     figure(3); hold on;
% % %     for yy = 1:length(ind3)
% % %         plot(singleCell(ind3(yy)).nf)
% % %     end
% % %     axis([0 60 1 7])
% % %     
% % %     figure(4); hold on;
% % %     for yy = 1:length(ind4)
% % %         plot(singleCell(ind4(yy)).nf)
% % %     end
% % %     axis([0 60 1 4])
    
% % % %% Function that finds the cell associated with trace and makes a movie of that cell
% % % img_dir = '\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20150710_PhosphateDeplet_ASOE_TRpanel\phosphateDeplet\';
% % % Nwells_singleCell = y_singlecell_cell(:,2);
% % % Nwells_t_singleCell = t_singlecell_cell(:,2);
% % % 
% % % for ii = 3:3%length(Nwells_singleCell{ii}) % loop through each well
% % %     
% % %     singleCell = Nwells_singleCell{ii};
% % %     tsingleCell = Nwells_t_singleCell{ii};
% % %     for kk = 3:3%length(singleCells)
% % %         
% % %   
% % %         img = imread(strcat(img_dir,'RFP_p7_t1.tiff'));
% % %         cX = singleCell(10).Cxloc;
% % %         cY = singleCell(10).Cyloc;
% % %         figure; imshow(img,[]); hold on; plot(cY(1),cX(1),'r.');
% % %         %for jj = 1:length(singleCell(kk).Cxloc)
% % %         %    img = imread(strcat(img_dir,filenames(jj).name));
% % %         %    cX = singleCell(kk).Cxloc(jj);
% % %         %    cY = singleCell(kk).Cyloc(jj);
% % %         %    figure(jj); imshow(img,[]); hold on; plot(cY,cX,'r.');
% % %         %end
% % %     end
% % % end

%% Plotting "Number of Peaks" Feature
TR_names1 = {'cad1', 'com2', 'crz1', 'dot6', 'maf1', 'msn2' ,'yap1', 'stb3', 'sko1', 'rtg3', 'pho4' ,'nrg2', 'msn4'};

numWells = length(out_Features_cell);
splot = ceil(sqrt(numWells));

% % create an empty map container
% peakNumObj = containers.Map;
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
    
%     % initialize the map container
%     for ff = 1:length(vals)
%         
%         display(strcat('This is ff: ', num2str(ff)))
%         
%         peaksAtValue = length(numPeaksVec(find(numPeaksVec == vals(ff))));
%         if isKey(peakNumObj, num2str(vals(ff))) % if the key already exists
%             display('1')
%             peakNumObj(num2str(vals(ff))) = [peakNumObj(num2str(vals(ff))), peaksAtValue];
% 
%         else % if key does not exist
%             display('2')
%             peakNumObj(num2str(vals(ff))) = peaksAtValue;
%         end
%     end

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

%% FIRST PEAK ONLY:
% SUBPLOT: MAX HEIGHT, TIME MAX HEIGHT, PULSE WIDTH, OFFSLOPE, ONSLOPE

TR_names1 = {'crz1', 'dot6', 'maf1', 'msn2' ,'yap1', 'stb3', 'sko1', 'rtg3', 'pho4' ,'nrg2', 'msn4'}; 

numWells = length(out_Features_cell);
%splot = ceil(sqrt(numWells));
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
    %subplot(splot,splot,jj)
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

%print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight.pdf','-dpdf')

maxHeightVecCell_1 = maxHeightVecCell;
maxHeightTimeVecCell_1 = maxHeightTimeVecCell;
pulseWidthVecCell_1 = pulseWidthVecCell;
offSlopeVecCell_1 = offSlopeVecCell;
onSlopeVecCell_1 = onSlopeVecCell;

%% SECOND PEAK ONLY:
% SUBPLOT: MAX HEIGHT, TIME MAX HEIGHT, PULSE WIDTH, OFFSLOPE, ONSLOPE
numWells = length(out_Features_cell);
%splot = ceil(sqrt(numWells));
pp =1;
%first = 2;
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
    %subplot(splot,splot,jj)
    numCells = length(out_Features_cell(jj).out_Features);
    for ii = 1:numCells

        % filter out the no peak data
        %if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
            
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
        %end

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

%% COMPARING FIRST AND SECOND PEAK FEATURES - plot over each other
figure(3); 
subplot(2,3,1); 
distributionPlot(maxHeightVecCell_1,'distWidth',0.9, 'color', 'g', 'addSpread',0, 'showMM',2,'xNames',TR_names1, 'yLabel', 'Maximum Height (nuc/cyt)'); title('Maximum Height'); hold on;
distributionPlot(maxHeightVecCell_2,'distWidth',0.9, 'color', 'c', 'addSpread',0, 'showMM',2);
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

%% COMPARISON OF FIRST AND SECOND PEAK FEATURES - ratios of 1st/2nd peak

% % %% THIRD PEAK FEATURES (sometimes there is no third peak)
% % % SUBPLOT: MAX HEIGHT, TIME MAX HEIGHT, PULSE WIDTH, OFFSLOPE, ONSLOPE
% % numWells = length(out_Features_cell);
% % %splot = ceil(sqrt(numWells));
% % pp =1;
% % %first = 2;
% % clear maxHeightVec
% % clear maxHeightVecCell
% % clear maxHeightTimeVec
% % clear maxHeightTimeVecCell
% % clear pulseWidthVec
% % clear pulseWidthVecCell
% % clear offSlopeVec
% % clear offSlopeVecCell
% % clear onSlopeVec
% % clear onSlopeVecCell
% % 
% % peakID = 3;
% % 
% % for jj = 3:numWells
% %     %subplot(splot,splot,jj)
% %     numCells = length(out_Features_cell(jj).out_Features);
% %     for ii = 1:numCells
% % 
% %         % filter out the no peak data
% %         %if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
% %             
% %             % only find the first peak
% %             if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == peakID %& first == 2
% % 
% %                 maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(peakID);
% %                 maxHeightTimeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(peakID);
% %                 pulseWidthVec(pp) = out_Features_cell(jj).out_Features(ii).cell.PulseWidth(peakID);
% %                 offSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OffSlope(peakID);
% %                 onSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OnSlope(peakID);
% %                 
% %                 %if ii >30 & ii <32;
% %                 %    figure;
% %                 %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
% %                 %    hold on;
% %                 %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
% %                 %end
% %                 
% %                 pp=pp+1;
% %             end
% %         %end
% % 
% %     end
% %     
% %     maxHeightVecCell{jj-2} = maxHeightVec;
% %     maxHeightTimeVecCell{jj-2} = maxHeightTimeVec;
% %     pulseWidthVecCell{jj-2} = pulseWidthVec;
% %     offSlopeVecCell{jj-2} = offSlopeVec;
% %     onSlopeVecCell{jj-2} = onSlopeVec;
% % end
% % 
% % 
% % % % %             % look at all the peaks
% % % % %             for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% % % % %                 maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(kk);
% % % % %                 pp=pp+1;
% % % % %             end
% %             
% % %print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight.eps','-depsc')
% % 
% % figure(4); 
% % subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',TR_names1, 'yLabel', 'Maximum Height (nuc/cyt)'); title('Maximum Height')
% % subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Maximum Height Time (min)'); title('Maximum Height Time');
% % subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Pulse Width (min)'); title('Pulse Width');
% % subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'Off Slope (nuc/cyt per 2min)'); title('Off Slope');
% % subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',TR_names1, 'yLabel', 'On Slope (nuc/cyt per 2min)'); title('On Slope');
% % suptitle('Third Peak Features');

% % %% Plotting histogram distributions of features
% % %clear
% % %%%%% base_dir = directory where single cell data is stored %%%%
% % %base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\processed_data_20150720\';
% % %%%%% specify the .mat file to load %%%%
% % %fname_load = '20150710_PhosphateDeplet_ASOE_TRpanel_20150720.mat';
% % %load([base_dir,fname_load],'out_Features_cell','TR_names', 'stress_type')
% % 
% % TR_names1 = {'cad1','com2', 'crz1', 'dot6', 'maf1', 'msn2' ,'yap1', 'stb3', 'sko1', 'rtg3', 'pho4' ,'nrg2', 'msn4'};
% % 
% % % Plotting: NumPeaks
% % figure(3);
% % numWells = length(out_Features_cell);
% % splot = ceil(sqrt(numWells));
% % clear numPeaksVec
% % clear numPeaksVecCell
% % for jj = 1:numWells
% %     subplot(splot,splot,jj)
% %     numCells = length(out_Features_cell(jj).out_Features)
% %     for ii = 1:numCells
% %         numPeaksVec(ii) = out_Features_cell(jj).out_Features(ii).cell.NumPeaks;
% %         
% % %         if ii >100 & ii < 112
% % %            figure
% % %            plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
% % %            hold on;
% % %            plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight', out_Features_cell(jj).out_Features(ii).cell.MaxHeight, 'ro');
% % %         end
% %         %if numPeaksVec(ii) == 1
% %         %   figure;
% %         %   plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
% %         %   hold on;
% %         %   plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight', out_Features_cell(jj).out_Features(ii).cell.MaxHeight, 'ro');
% %         %end
% %     end
% %     hist(numPeaksVec, sqrt(length(numPeaksVec)));
% %     %%%% hard coded range %%%%
% %     axis([0 3 0 200]);
% %     title(TR_names1{jj})
% %     
% %     numPeaksVecCell{jj} = numPeaksVec;
% % end
% % 
% % figure; distributionPlot(numPeaksVecCell,'distWidth',0.9,'addSpread',1,'showMM',5,'xNames',TR_names1, 'yLabel', 'Number Of Peaks')
% % title('Number Peaks: pH Stress');
% % %rotateXlabels(handles,90)
% % print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\NumPeaks.pdf','-dpdf')

%%%%%%%%%It does looks like there are more than 1 peak in a particular
%%%%%%%%%trace, but no more than 3 peaks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%
% % % Plotting: MaxHeight (1st peak OR all peaks)
% % figure(4);
% % numWells = length(out_Features_cell);
% % splot = ceil(sqrt(numWells));
% % pp =1;
% % first = 0;
% % clear maxHeightVec
% % clear maxHeightVecCell
% % for jj = 1:numWells
% %     subplot(splot,splot,jj)
% %     numCells = length(out_Features_cell(jj).out_Features);
% %     for ii = 1:numCells
% %         %if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == 0
% %             %maxHeightVec(pp) = 0;
% %             %pp= pp+1;
% %         %else
% %         % filter out the no peak data
% %         if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
% %             % look at all the peaks
% %             if first == 0;
% %                 for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% %                     maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(kk);
% %                     pp=pp+1;
% %                 end
% %             % look at only the first peak
% %             else
% %                 maxHeightVec(pp) = out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1);
% %                 
% %                 %if ii >30 & ii <32;
% %                 %    figure;
% %                 %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
% %                 %    hold on;
% %                 %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
% %                 %end
% %                 
% %                 pp=pp+1;
% %             end
% %         end
% %         %end
% %     end
% %     hist(maxHeightVec, sqrt(length(maxHeightVec)));
% %     %%%% hard coded range %%%%
% %     axis([0 12 0 200]);
% %     title(TR_names{jj})
% %     
% %     maxHeightVecCell{jj} = maxHeightVec;
% % end
% % 
% % %print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight.eps','-depsc')
% % 
% % figure; distributionPlot(maxHeightVecCell,'distWidth',0.9,'addSpread',1,'showMM',5,'xNames',TR_names1, 'yLabel', 'Maximum Height (nuc/cyt)');
% % title('Maximum Height: pH Stress');
% % print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight.pdf','-dpdf')
% % %%%%%%%%%% most of the max heights are similar, but some have higher tails
% % %%%%%%%%%% %%%%%%%%%%%%%%%%%
% % %%
% % % Plotting: MaxHeightTime
% % figure(5);
% % numWells = length(out_Features_cell);
% % splot = ceil(sqrt(numWells));
% % pp =1;
% % first = 0;
% % clear maxHeightTimeVec
% % clear maxHeightTimeVecCell
% % for jj = 1:numWells
% %     subplot(splot,splot,jj)
% %     numCells = length(out_Features_cell(jj).out_Features);
% %     for ii = 1:numCells
% %         %if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == 0
% %         %    maxHeightTimeVec(pp) = 0;
% %         %    pp= pp+1;
% %         %else
% %         if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
% %             if first == 0; % Height Time for multiple peaks in the same trace
% %                 for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% %                     maxHeightTimeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(kk);
% %                     pp=pp+1;
% %                 end
% %             else
% %                 maxHeightTimeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1);
% %                 
% % %                 if maxHeightTimeVec(pp) > 20 & maxHeightTimeVec(pp) < 40
% % %                    figure;
% % %                    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
% % %                    hold on;
% % %                    plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
% % %                 elseif maxHeightTimeVec(pp) > 50 & maxHeightTimeVec(pp) < 70
% % %                    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace); 
% % %                    hold on;
% % %                    plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
% % %                 elseif maxHeightTimeVec(pp) > 80 & maxHeightTimeVec(pp) < 100
% % %                    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes', out_Features_cell(jj).out_Features(ii).cell.TimeTrace);
% % %                    hold on;
% % %                    plot(out_Features_cell(jj).out_Features(ii).cell.TimeMaxHeight(1)', out_Features_cell(jj).out_Features(ii).cell.MaxHeight(1), 'ro');
% % %                 end
% %                 
% %                 pp=pp+1;
% %             end
% %         end
% %         %end
% %     end
% %     hist(maxHeightTimeVec, sqrt(length(maxHeightTimeVec)));
% %     %%%% hard coded range %%%%
% %     axis([0 140 0 100]);
% %     title(TR_names{jj})
% %     
% %     maxHeightTimeVecCell{jj} = maxHeightTimeVec;
% % end
% % 
% % %print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeightTime.eps','-depsc')
% % figure; distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',1,'showMM',5,'xNames',TR_names1, 'yLabel', 'Max Height Time (min)');
% % title('Maximum Height Time: pH Stress');
% % print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeightTime.pdf','-dpdf')
%%%%%%%%%%%%%%% the times of max height actually differ based on the TR,
%%%%%%%%%%%%%%% and for some TRs, there are multiple times of peaking %%%%%
% % %%
% % % Plotting: PulseWidth
% % figure(6);
% % numWells = length(out_Features_cell);
% % splot = ceil(sqrt(numWells));
% % pp =1;
% % first = 0;
% % clear pulseWidthVec
% % clear pulseWidthVecCell
% % for jj = 1:numWells
% %     subplot(splot,splot,jj)
% %     numCells = length(out_Features_cell(jj).out_Features);
% %     for ii = 1:numCells
% %         %if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == 0
% %         %    pulseWidthVec(pp) = 0;
% %         %    pp= pp+1;
% %         %else
% %         if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
% %             if first == 0; % Height Time for multiple peaks in the same trace
% %                 for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% %                     pulseWidthVec(pp) = out_Features_cell(jj).out_Features(ii).cell.PulseWidth(kk);
% %                     pp=pp+1;
% %                 end
% %             else
% %                 pulseWidthVec(pp) = out_Features_cell(jj).out_Features(ii).cell.PulseWidth(1);
% %                 
% %                 %if ii > 0 & ii < 11
% %                 %    figure;
% %                 %    %plot(out_Features_cell(jj).out_Features(ii).cell.TimeTrace); hold on;
% %                 %    plot(out_Features_cell(jj).out_Features(ii).cell.TimeTraceTimes, out_Features_cell(jj).out_Features(ii).cell.TimeTrace); hold on; 
% %                 %    plot(out_Features_cell(jj).out_Features(ii).cell.OffSlopeTime(1),2,'r.');
% %                 %    plot(out_Features_cell(jj).out_Features(ii).cell.OnSlopeTime(1), 2, 'r.');
% %                 %end
% %                 
% %                 pp=pp+1;
% %             end
% %         end
% %     end
% %     hist(pulseWidthVec, sqrt(length(pulseWidthVec)));
% %     %%%% hard coded range %%%%
% %     axis([0 120 0 100]);
% %     title(TR_names{jj})
% %     
% %     pulseWidthVecCell{jj} = pulseWidthVec;
% % end
% % 
% % %print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\PulseWidth.eps','-depsc')
% % 
% % figure; distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',1,'showMM',5,'xNames',TR_names1, 'yLabel', 'Pulse Width (min)');
% % title('Pulse Width: pH Stress');
% % print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\PulseWidth.pdf','-dpdf')

%%%%%%% the pulse width does not look like it is helpful because of how
%%%%%%% spread it is %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%
% % % Plotting: OnSlope
% % figure(7);
% % numWells = length(out_Features_cell);
% % splot = ceil(sqrt(numWells));
% % pp =1;
% % first = 0;
% % clear onSlopeVec
% % clear onSlopeVecCell
% % for jj = 1:numWells
% %     subplot(splot,splot,jj)
% %     numCells = length(out_Features_cell(jj).out_Features);
% %     for ii = 1:numCells
% %         %if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == 0
% %         %    onSlopeVec(pp) = 0;
% %         %    pp= pp+1;
% %         %else
% %         if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
% %             if first == 0; % Height Time for multiple peaks in the same trace
% %                 for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% %                     onSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OnSlope(kk);
% %                     pp=pp+1;
% %                 end
% %             else
% %                 onSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OnSlope(1);
% %                 pp=pp+1;
% %             end
% %         end
% %     end
% %     hist(onSlopeVec, sqrt(length(onSlopeVec)));
% %     %%%% hard coded range %%%%
% %     axis([0 1 0 250]);
% %     title(TR_names{jj})
% %     
% %     onSlopeVecCell{jj} = onSlopeVec;
% % end
% % 
% % %print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\OnSlope.eps','-depsc')
% % figure; distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',1,'showMM',5,'xNames',TR_names1, 'yLabel', 'On Slope (d(nuc/cyt)/dt) dt = every 2 min');
% % title('On Slope: pH Stress');
% % 
% % print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\OnSlope.pdf','-dpdf')
%%%%%%%% much of the on slope has a single peak with varying tail %%%%%%%%%
% % %%
% % % Plotting: OffSlope
% % figure(8);
% % numWells = length(out_Features_cell);
% % splot = ceil(sqrt(numWells));
% % pp =1;
% % first = 0; % this is important - tells you whether a single trace has a multimodal shape
% % clear offSlopeVec
% % clear offSlopeVecCell
% % for jj = 1:numWells
% %     subplot(splot,splot,jj)
% %     numCells = length(out_Features_cell(jj).out_Features);
% %     for ii = 1:numCells
% %         %if out_Features_cell(jj).out_Features(ii).cell.NumPeaks == 0
% %         %    offSlopeVec(pp) = 0;
% %         %    pp= pp+1;
% %         %else
% %         if out_Features_cell(jj).out_Features(ii).cell.NumPeaks ~= 0
% %             if first == 0; % Height Time for multiple peaks in the same trace
% %                 for kk = 1:out_Features_cell(jj).out_Features(ii).cell.NumPeaks
% %                     offSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OffSlope(kk);
% %                     pp=pp+1;
% %                 end
% %             else
% %                 offSlopeVec(pp) = out_Features_cell(jj).out_Features(ii).cell.OffSlope(1);
% %                 pp=pp+1;
% %             end
% %         end
% %     end
% %     hist(offSlopeVec, sqrt(length(offSlopeVec)));
% %     %%%% hard coded range %%%%
% %     axis([-0.7 0 0 200]);
% %     title(TR_names{jj})
% %     
% %     offSlopeVecCell{jj} = offSlopeVec;
% %     
% % end
% % 
% % %print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\OffSlope.eps','-depsc')
% % figure; distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',1,'showMM',5,'xNames',TR_names1, 'yLabel', 'Off Slope (d(nuc/cyt)/dt) dt = every 2 min');
% % title('Off Slope: pH Stress')
% % print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\OffSlope.pdf','-dpdf')

%%%%%%%%%%%%%%%% looks like they all have the same off slope. kind of similar to on slope %%%%%%%%%%%%%%%%%

% The fact that these have such a wide spread on some of the features makes
% it hard to be useful

% When features are correlated, then it means that could be regulated by
% the same process, AND that you only have to consider one of the features.

% WRITE ANOTHER SCRIPT TO PLOT TRS VERSUS STRESS INPUTS
%%%%
% NORMALIZE THE DATA
% CLUSTER TO EXTRACT THE FEATURES IN AN UNSUPERVISED WAY
% HOLT-WINTERS ALGORITHM

%%%%
% NumPeak bar graph
%% Plotting correlations among features

% Max Height VS Time Max Height
figure(9)
for ii = 1:length(maxHeightTimeVecCell)
    ii
    currMaxHeight = maxHeightVecCell{ii};
    currMaxHeightTime = maxHeightTimeVecCell{ii};
    subplot(ceil(sqrt(length(maxHeightTimeVecCell))), ceil(sqrt(length(maxHeightTimeVecCell))),ii)
    plot(currMaxHeight',currMaxHeightTime,'.');
    %refline(1,0)
    title(TR_names{ii})
    axis([1.5 10 0 150])
end

print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight_vs_MaxHeightTime.eps','-depsc')
% Max Height VS Pulse Width
figure(10)
for ii = 1:length(maxHeightVecCell)
    ii
    currMaxHeight = maxHeightVecCell{ii};
    currPulseWidth = pulseWidthVecCell{ii};
    subplot(ceil(sqrt(length(maxHeightVecCell))), ceil(sqrt(length(maxHeightVecCell))),ii)
    plot(currMaxHeight,currPulseWidth,'.');
    title(TR_names{ii})
end
print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight_vs_PulseWidth.eps','-depsc')
% Max Height VS On Slope
figure(11)
for ii = 1:length(maxHeightVecCell)
    ii
    currMaxHeight = maxHeightVecCell{ii};
    currOnSlope = onSlopeVecCell{ii};
    subplot(ceil(sqrt(length(maxHeightVecCell))), ceil(sqrt(length(maxHeightVecCell))),ii)
    plot(currMaxHeight,currOnSlope,'.');
    title(TR_names{ii})
end
print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\MaxHeight_vs_OnSlope.eps','-depsc')
% Max Height VS Off Slope
figure(12)
for ii = 1:length(maxHeightVecCell)
    ii
    currMaxHeight = maxHeightVecCell{ii};
    currOffSlope = offSlopeVecCell{ii};
    subplot(ceil(sqrt(length(maxHeightVecCell))), ceil(sqrt(length(maxHeightVecCell))),ii)
    plot(currMaxHeight,currOffSlope,'.');
    title(TR_names{ii})
end
print('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150710_PhosphateDeplet_ASOE_TRpanel\plots_20150720\OffSlope.eps','-depsc')

