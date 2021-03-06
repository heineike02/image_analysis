% SCRIPT to generate test plots - TEST PLOTS

% run the find cells script for only 1 well -- this is done with
% SCY_img_analysis_20150820_Dot6_7StressInputs_1SDC.m for example

addpath('/Users/susanychen/GITREPOS/image_analysis2');

% %plot the tracks on top of image - make sure reasonable cells were found
im = imread('RFP_p1_t32.tiff');
figure; imshow(im,[]);
hold on;
for i = 1:length(all_tracks_vec{1}.one)
plot(all_tracks_vec{1}.one(i).Cyloc(1),all_tracks_vec{1}.one(i).Cxloc(1), 'g.'); hold on;
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
fname_load = 'NLSvar11_doseresp_freq_20150922';
val_to_plot = 'nf';

% initialize
base_dir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20150922_NLSvar11_doseresp_freq\';
experi_name = fname_load;
load([base_dir,experi_name])

phases = {'one','two','three','four','five','six', 'seven','eight','nine','ten','eleven'};

well_order = {'No Light','25auBL 3min','50auBL 3min','100auBL 3min','150auBL 3min','255auBL 3min','10minPeriod@255','6minPeriod@255','4minPeriod@255','2minPeriod@255','1minPeriod@255'};
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

defineLight(1).light = [2*ones(1,8),2*ones(1,12), 2*ones(1,12)];
defineLight(2).light = [2*ones(1,8),(2+(25/255))*ones(1,12), 2*ones(1,12)];
defineLight(3).light = [2*ones(1,8),(2+(50/255))*ones(1,12), 2*ones(1,12)];
defineLight(4).light = [2*ones(1,8),(2+(100/255))*ones(1,12), 2*ones(1,12)];
defineLight(5).light = [2*ones(1,8),(2+(150/255))*ones(1,12), 2*ones(1,12)];
defineLight(6).light = [2*ones(1,8),(2+(255/255))*ones(1,12), 2*ones(1,12)];
defineLight(7).light = [3*ones(1,10), 2*ones(1,10), 3*ones(1,10), 2*ones(1,10)];
defineLight(8).light = [3*ones(1,6), 2*ones(1,6), 3*ones(1,6), 2*ones(1,6)];
defineLight(9).light = [3*ones(1,8), 2*ones(1,8), 3*ones(1,8), 2*ones(1,8)];
defineLight(10).light = [3*ones(1,4), 2*ones(1,4), 3*ones(1,4), 2*ones(1,4)];
defineLight(11).light = [3*ones(1,3), 2*ones(1,3), 3*ones(1,3), 2*ones(1,3)];

%defineLight(1).light = [3*ones(1,6),4*ones(1,6), 3*ones(1,6), 4*ones(1,6), 3*ones(1,6)];

theAxis{1} = [0 8 1.5 3];
theAxis{7} = [0 20 1.5 3];
theAxis{8} = [0 12 1.5 3];
theAxis{9} = [0 8 1.5 3];
theAxis{10} = [0 4 1.5 3];
theAxis{11} = [0 2 1.5 3];
for jj = 1:Nfolders
    jj
    all_tracks = all_tracks_vec{1}.(phases{jj});
    all_times = all_times_vec{1}.(phases{jj});

    %subplot(num_subplots,num_subplots,jj);
    subplot(3,4,jj);
    
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
    text(1,2.5,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    if jj <=6
        axis(theAxis{1});
    else
        axis(theAxis{jj});
    end
    
    title(sample_order{jj})
end
suptitle(save_img_file(1).filename);
%print([base_dir, save_img_file(1).filename, ' 20150922 NLSvar11 doseresp freq','.eps'], '-depsc')

%% plotting mean and standard deviation
%% %% Plot the mean with standard deviation
figure(1)
clf

hold on
%%%% specify channel name %%%%
channel = 'RFP'; %'RFP'
cmap_RFP = [ 0,0,1;  %Blue
    0,1,0; 
    1 0 0; 
    0,1,1;  
    1,0,1; 
    1 1 0;
    0,1,0; % start of next set 
    0,1,0;
    1,0,0;%red
    1 0 1; %color
    1 0 1;
    1 0 1;
    ];  
plot_params = {'linewidth',1.5,'LineStyle','-'};

setAxisValues{7} = [0 20 1.5 2.5];
setAxisValues{8} = [0 12 1.5 2.5];
setAxisValues{9} = [0 8 1.5 2.5];
setAxisValues{10} = [0 4 1.5 2.5];
setAxisValues{11} = [0 2 1.5 2.5];

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
% for plotting specific subplots

% extract the max/min difference for a frequency vs localization plot (20170215)

%meanHighLow = [];
%stdHighLow = [];
plotNums = [8,9,10];

%valsStruct = struct();
for jj = 1:length(plotNums)%7:11%Nfolders
    jj
    % all_tracks = all_tracks_vec{1}.(phases{jj});
    % all_times = all_times_vec{1}.(phases{jj});
    all_tracks = all_tracks_vec{1}.(phases{plotNums(jj)});
    all_times = all_times_vec{1}.(phases{plotNums(jj)});

    %subplot(num_subplots,num_subplots,jj);
    %subplot(1,5,jj-6);
    figure(plotNums(jj))
    
    % loop through single cells
    %hold all
    xlabel('time(min)')
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
    p= plot_meanvalues_syc(all_times,all_tracks,channel,color_val,0,'nf', setAxisValues{jj}, 'plot_params',plot_params);
    %[fig_out, timevals, mean_val, std_val] = plot_meanvalues_syc(all_times,all_tracks,channel,color_val,0,'nf', setAxisValues{jj}, 'plot_params',plot_params);
    
    
    %[maxVal, indxMax] = max(mean_val); [minVal, indxMin]=min(mean_val(postTime(jj-6):end));
    %meanHighLow(jj) = maxVal - minVal;
    %stdHighLow(jj) = sqrt(std_val(indxMax).^2 + std_val(postTime(jj-6)+indxMin).^2);

    %plot(all_times, defineLight(1).light, 'c-.', 'linewidth',4);
    
    %valsStruct(jj-6).mean = mean_val;
    %valsStruct(jj-6).std = std_val;
    
    grid off
    %totCellNum = length(all_tracks);
    %text(0.5,2.25,['CellNum: ',num2str(totCellNum)])
    
    %%%% specify x and y range of plot %%%%
    %if jj <=6
        %axis(theAxis{1});
        %xlim([0 12])
    %else
        %axis(theAxis{jj});
    %end
    
    title(sample_order{plotNums(jj)})
    %title(sample_order{jj})
end
suptitle(save_img_file(1).filename);

% %% extract the max/min difference for a frequency vs localization plot (20170215)
% % imaging frequency is every 30secs
% phases = {'seven','eight','nine','ten','eleven'};
% begTime2Sample = [5,5,4,2,1]; %indices
% endTime2Sample = [15,10,7,4,2]; %indices
% 
% highLowDiffMed = [];
% highLowDiffStd = [];
% for jj = 1:5
%     %jj
%     all_tracks = all_tracks_vec{1}.(phases{jj});
%     all_times = all_times_vec{1}.(phases{jj});
%     
%     %highCount = []; % starts with light on, so localization is high
%     %lowCount = [];
%     %highLowDiff = [];
%     rangeDiff = [];
%     for ii = 1:length(all_tracks)
%         %ii
%         tracktime1 = all_tracks(ii).times;
%         track1 = all_tracks(ii).nf.RFP;
%         
%         maxSingle = max(track1); minSingle = min(track1);
%         rangeDiff(ii) = maxSingle - minSingle;
%         %if find(tracktime1 == begTime2Sample(jj)) & find(tracktime1 == endTime2Sample(jj))
%             %figure; plot(track1.RFP)
%             %indxHigh = find(tracktime1 == begTime2Sample(jj));
%             %indxLow = find(tracktime1 == endTime2Sample(jj));
%             %highCount(ii) = track1.RFP(indxHigh);
%             %lowCount(ii) = track1.RFP(indxLow);
%             %highLowDiff(ii) = track1.RFP(indxHigh)-track1.RFP(indxLow);
%         %end
%     end
%     %highLowDiffMed(jj) = nanmedian(highLowDiff);
%     %highLowDiffStd(jj) = nanstd(highLowDiff);
%     highLowDiffMed(jj) = nanmedian(rangeDiff);
%     highLowDiffStd(jj) = nanstd(rangeDiff);
% end

%% extracting max min for freq plot
%7-11
meanVals = [0.45,0.55,0.48,0.43,0.17];
maxIndx = [3 3 5 5 5];
minIndx = [16 12 17 9 7];
stdVals = [];
for i = 1:5
    stdVals(i) = sqrt(valsStruct(i).std(minIndx(i)).^2+valsStruct(i).std(maxIndx(i)).^2);
end
figure; shadedErrorBar([1/10, 1/6, 1/4, 1/2, 1/1],meanVals, stdVals); %errorbar([1/10, 1/6, 1/4, 1/2, 1/1],meanVals, stdVals);
axis([0 1 0 0.8]);