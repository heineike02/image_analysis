function [] = SYC_filter_extraction_plotting_SC_template_2(base_dir, experi_name, well_order, phases, save_img_file,plotVars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Script plots the same summary plots as for each of the experiments

% Example:
%   base_dir = '/Users/susanychen/Downloads/TempMatlabData20150808/';
%   experi_name = '20150629_GD_doseresp_Msn2Msn4Maf1Stb3_20150720_filter';
%   phases = {'Pre','Post'};
%   save_img_file = structure of filenames for each of the plots generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
base_dir = base_dir;
experi_name = experi_name;
load([base_dir,experi_name])
%TR_names_long(9:12);
%TR_names_short(9:12);
timeVecs = all_times_filt_vec;
trackVecs = all_tracks_filt_vec;
featStructCell = out_Features_cell; % just access the peak vs no peak field
phases = phases;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Mean Time Traces - filtered - Msn2
% figure(1)
% clf
% %hold on
% clear mean_val_cell; clear std_val_cell;
% 
% condNames = {'GD', 'GD1', 'OS', 'ZM', 'PH'};
% %%%% specify the subfolders present in the experiment %%%%
% phases = {'Pre', 'Post'};
% %%%% leave this blank if 1 channel, write channel name if more than 1 %%%%
% channel = '';
% %%%% specify color of lines in plot %%%%
% color_val = 'b'; %cmap(jj,:);
% %%%% specify plot line properties %%%%
% plot_params = {'linewidth',1.5,'LineStyle','-'};
% %%%% 'nf' = nuclear localization; 'nmi' = average fluorescence %%%%
% val_to_plot = 'nf';
% 
% Nwells = length(Msn2_5cond_time);
% 
% %num_subplots = ceil(sqrt(Nwells));
% 
% % "all_tracks_vec" and "all_times_vec" are both 1xNwells cells
% for theWells = 1:Nwells % loop through each well
%     
%     subplot(2,3,theWells);
%             
%     all_tracks = Msn2_5cond_track{theWells};
%     all_times = Msn2_5cond_time{theWells};
%     for pH = 1:length(phases) % loop through pre and post phases
%         tracks = all_tracks.(phases{pH}); length(tracks)
%         timevals = all_times.(phases{pH});
% 
%         [p, timevals_cell{theWells,pH},mean_val_cell{theWells,pH}, std_val_cell{theWells,pH}] = plot_meanvalues_syc(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params); % either line plot or line plot with errorbars
%         hold on;
%         title(condNames{theWells})
%         %%%% specify the x and y range of the plot %%%%
%         ylim([1 5])
%         xlim([0 100])
%         %axis([0 130 1 6])
%     end
%     
% end
% suptitle('Msn2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SC.MSN2-RFP plot individual cells (filtered, not smoothed)
% - variance of the time trace (captures presence of responders and non and
% different types of profiles)
% Msn2
figure(1)
clf
clear t_singlecell_cell; clear y_singlecell_cell;
hold on
%%%% specify channel name %%%%
channel = 'RFP'; %'RFP'

well_order = well_order;
%condNames = {'GD', 'GD1', 'OS', 'ZM', 'PH'};
 
Nwells = length(trackVecs);
num_subplots = ceil(sqrt(Nwells));

for jj = 1:Nwells
    jj
    all_tracks = trackVecs{jj};
    all_times = timeVecs{jj};
    %subplot(1,2,jj)
    subplot(num_subplots,num_subplots,jj);
    for ph = 1: length(phases) % loop through pre and post % is this defined??
        
        tracks = all_tracks.(phases{ph}); size(tracks)
        timevals = all_times.(phases{ph}); %size(timevals)
        [p] = plot_individual_values(timevals,tracks,channel,'nf');
        hold on;
        
        t_singlecell_cell{jj,ph} = timevals;
        y_singlecell_cell{jj,ph} = tracks;
    end
    
    totCellNum = featStructCell(jj).hasPeakNum+featStructCell(jj).zeroPeakNum;
    text(1,1,['CellNum: ',num2str(totCellNum)])
    %%%% specify x and y range of plot %%%%
    axis(plotVars.axis)
    title(well_order{jj})
end
suptitle(save_img_file(1).filename);
%print([base_dir, save_img_file(1).filename,'.eps'], '-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotting single cell plots with peaks indicated
% 
% % Msn2 environmental conditions
% for bb = 1:length(Msn2_5cond_feat) % this is each environmental condition
%     currEnvCond = Msn2_5cond_feat(bb).Msn2;
%     figure;
%     if length(currEnvCond) > 81 
%         for aa = 1:81 % go through every single cell
%             subplot(9,9,aa);
%             currCell = currEnvCond(aa).cell;
%             plot(currCell.TimeTraceTimes, currCell.TimeTrace); 
%             if currCell.NumPeaks ~= 0
%                 hold on;
%                 plot(currCell.TimeMaxHeight, currCell.MaxHeight, 'r.');
%                 if ~isempty(currCell.OnSlopeTime)
%                     plot([currCell.OnSlopeTime], [2], 'g.');
%                 end
%                 if ~isempty(currCell.OffSlopeTime)
%                     plot([currCell.OffSlopeTime], [2], 'g.');
%                 end
%             end
%             axis([0 100 0 5])
%         end
%     else 
%         numSub = ceil(sqrt(length(currEnvCond)));
%         for aa = 1:length(currEnvCond);
%             subplot(numSub,numSub,aa);
%             currCell = currEnvCond(aa).cell;
%             plot(currCell.TimeTraceTimes, currCell.TimeTrace); 
%             if currCell.NumPeaks ~= 0
%                 hold on;
%                 plot(currCell.TimeMaxHeight, currCell.MaxHeight, 'r.');
%                 if ~isempty(currCell.OnSlopeTime)
%                     plot([currCell.OnSlopeTime], [2], 'g.');
%                 end
%                 if ~isempty(currCell.OffSlopeTime)
%                     plot([currCell.OffSlopeTime], [2], 'g.');
%                 end
%             end
%             axis([0 100 0 5])
%         end
%     end
%     suptitle(['Msn2:', condNames1{bb}])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting "Number of Peaks" Feature
%condNames1 = {'GD', 'GD1', 'OS', 'ZM', 'PH'};

numWells = length(featStructCell);
splot = ceil(sqrt(numWells));

clear uniqueVals
clear peakNumVecCell

for jj = 1:numWells
    
    display(strcat('This is jj: ', num2str(jj)))
    
    numCells = length(featStructCell(jj).out_Features);
    
    for ii = 1:numCells
        numPeaksVec(ii) = featStructCell(jj).out_Features(ii).cell.NumPeaks;
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
set(gca, 'xticklabel', well_order)
%xlabel('TRs'); ylabel('Number of Cells (count)')
legend(strsplit(num2str(diffPeaks)))
%title('Proportion of Number of Peaks Msn2')
%print([base_dir, save_img_file(2).filename], '-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST PEAK ONLY:
% SUBPLOT: MAX HEIGHT, TIME MAX HEIGHT, PULSE WIDTH, OFFSLOPE, ONSLOPE

%condNames = {'GD', 'GD1', 'OS', 'ZM' ,'PH'};

numWells = length(trackVecs);

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
    
    numCells = length(featStructCell(jj).out_Features);
    for ii = 1:numCells
        %ii
        % filter out the no peak data
        if featStructCell(jj).out_Features(ii).cell.NumPeaks ~= 0
           
            % only find the first peak (could have more than 1 peak)!
            if featStructCell(jj).out_Features(ii).cell.NumPeaks >= 1 & ~isempty(featStructCell(jj).out_Features(ii).cell.MaxHeight(1))
                currMsn2 = featStructCell(jj).out_Features;
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
subplot(2,3,1); distributionPlot(maxHeightVecCell,'distWidth',0.9, 'addSpread',0, 'showMM',5,'xNames',well_order, 'yLabel', ' (nuc/cyt)'); title('Maximum Height')
subplot(2,3,2); distributionPlot(maxHeightTimeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',well_order, 'yLabel', '(min)'); title('Maximum Height Time');
subplot(2,3,3); distributionPlot(pulseWidthVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',well_order, 'yLabel', '(min)'); title('Pulse Width');
subplot(2,3,4); distributionPlot(offSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',well_order, 'yLabel', '(nuc/cyt per 2min)'); title('Off Slope');
subplot(2,3,5); distributionPlot(onSlopeVecCell,'distWidth',0.9,'addSpread',0,'showMM',5,'xNames',well_order, 'yLabel', '(nuc/cyt per 2min)'); title('On Slope');
suptitle(save_img_file(3).filename);

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 10 10])
%print([base_dir, save_img_file(3).filename, '.eps'], '-depsc')

end
