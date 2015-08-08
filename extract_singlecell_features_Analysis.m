function Features = extract_singlecell_features_Analysis(tsingleCells, singleCells, var_to_plot, smoothFn, peakFindr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extracts features of a single cell trace: You can choose
% which features to output
% - number of peaks
% - max height(s)
% - time(s) of max height
% - rise time (on slope)
% - duration of pulse
%
% One Important Caveat is that this function quantifies transient PULSE
% characteristics. It does not deal well with FLAT or STEP time series
% data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% isolate the values for the analysis
singleCellTrace = singleCells.(var_to_plot);
% isolate the times for analysis
singleCellTraceTime = singleCells.times;
singleCellTraceTime = tsingleCells(singleCellTraceTime);

% smoothed single cell trace
if smoothFn.flag == 1
    %figure; 
    %plot(singleCellTraceTime, singleCellTrace);
    singleCellTrace = smooth(singleCellTrace, smoothFn.span, smoothFn.method);
    %hold on; plot(singleCellTraceTime, singleCellTrace,'r');
end

% Features: Number of peaks
% Method: Use the pre-packaged peakFinder function
sel = peakFindr.Sel;
thresh = peakFindr.Thresh;
extrema = 1;
incl_endpt = 0;
stdMinus = peakFindr.StdMinus;
[peakLoc, peakMag] = peakfinder_20150730(singleCellTrace, sel, thresh, extrema, incl_endpt, stdMinus);

[peakLoc, indc]=unique(peakLoc);
peakMag = peakMag(indc);

% test
%figure(1); plot(singleCellTraceTime, singleCellTrace); hold on;
%plot(singleCellTraceTime(peakLoc), peakMag, 'r.');

% take out peaks that are at the beginning and the ends
if peakLoc == 1
    peakLoc(1) = [];
    peakMag(1) = [];
end
if peakLoc == length(singleCellTrace)
    peakLoc(end) = [];
    peakMag(end) = [];
end

Features.NumPeaks = length(peakLoc);

% Features: Max Height and Time Max Height
% Method: Use the pre-packaged peakFinder function
Features.MaxHeight = peakMag;
Features.TimeMaxHeight = singleCellTraceTime(peakLoc);

% Features: time trace and the times
Features.TimeTrace = singleCellTrace;
Features.TimeTraceTimes = singleCellTraceTime;

% Features: the ON slope of a pulse
% Method: For each peak, look at the maximum slope to the left of it
% Features: duration of pulse
% Method: Find difference between times of OnSlope and OffSlope
numPeaks = length(peakLoc);
if numPeaks == 1

    onSlope = diff(singleCellTrace(1:peakLoc));
    [maxOnSlope, maxIndc]=max(onSlope);

    
    
    offSlope = diff(singleCellTrace(peakLoc:end));
    [minOffSlope, minIndc]=min(offSlope);
    
    % if there is no OffSlope or PulseWidth associated with peak,
    % then put NaN. Disregard peak.
    if isempty(minOffSlope) | isempty(maxOnSlope)
        Features.OffSlope = NaN;
        Features.OffSlopeTime = NaN;
        
        Features.OnSlope = NaN;
        Features.OnSlopeTime = NaN;

        Features.PulseWidth = NaN;
    else
%         % test
%         figure(2); plot(singleCellTraceTime, singleCellTrace); hold on; 
%         plot(singleCellTraceTime(maxIndc),singleCellTrace(maxIndc),'r.');   
%         plot(singleCellTraceTime(peakLoc+minIndc),singleCellTrace(peakLoc+minIndc),'r.');
            
        Features.OffSlope = minOffSlope;
        Features.OffSlopeTime = singleCellTraceTime(peakLoc+minIndc);
        Features.OnSlope = maxOnSlope;
        Features.OnSlopeTime = singleCellTraceTime(maxIndc);
        Features.PulseWidth = singleCellTraceTime(peakLoc+minIndc) - singleCellTraceTime(maxIndc);
    end
    
elseif numPeaks > 1 % if there are more than 1 peak

    for ii = 1:numPeaks
        if ii == 1 % first peak

            onSlope = diff(singleCellTrace(1:peakLoc(ii)));
            [maxOnSlope, maxIndc]=max(onSlope);
            
            
            offSlope = diff(singleCellTrace(peakLoc(ii):peakLoc(ii+1)));
            [minOffSlope, minIndc]=min(offSlope);
  
            % if there is no OffSlope or PulseWidth associated with peak,
            % then put NaN. Disregard peak.
            if isempty(minOffSlope)|isempty(maxOnSlope)
                temp_OffSlope(ii) = NaN;
                temp_OffSlopeTime(ii) = NaN;
                
                temp_OnSlope(ii) = NaN;
                temp_OnSlopeTime(ii) = NaN;
                
                temp_PulseWidth(ii) = NaN;
            else
%                 % test
%                 figure(2); plot(singleCellTraceTime, singleCellTrace); hold on; 
%                 plot(singleCellTraceTime(maxIndc),singleCellTrace(maxIndc),'r.');
%                 plot(singleCellTraceTime(peakLoc(ii)+minIndc),singleCellTrace(peakLoc(ii)+minIndc),'r.');
                
                temp_OffSlope(ii) = minOffSlope;
                temp_OffSlopeTime(ii) = singleCellTraceTime(peakLoc(ii)+minIndc);
                temp_OnSlope(ii) = maxOnSlope;
                temp_OnSlopeTime(ii) = singleCellTraceTime(maxIndc);
                temp_PulseWidth(ii) = singleCellTraceTime(peakLoc(ii)+minIndc) - singleCellTraceTime(maxIndc);
            end
            
        end
        
        if ii+1 > numPeaks % last peak - note that peakLoc = index. 
   
            onSlope = diff(singleCellTrace(peakLoc(ii-1):peakLoc(ii)));
            [maxOnSlope, maxIndc]=max(onSlope);
            

            temp_OnSlope(ii) = maxOnSlope;
            temp_OnSlopeTime(ii) = singleCellTraceTime(peakLoc(ii-1)+maxIndc);
            
            offSlope = diff(singleCellTrace(peakLoc(ii):end));
            [minOffSlope, minIndc]=min(offSlope);
            

            % if there is no OffSlope or PulseWidth associated with peak,
            % then put NaN. Disregard peak.
            if isempty(minOffSlope)|isempty(maxOnSlope)
                temp_OffSlope(ii) = NaN;
                temp_OffSlopeTime(ii) = NaN;
                temp_OnSlope(ii) = NaN;
                temp_OnSlopeTime(ii) = NaN;
                temp_PulseWidth(ii) = NaN;
            else
%                 % test
%                 figure(2); plot(singleCellTraceTime, singleCellTrace); hold on; 
%                 plot(singleCellTraceTime(peakLoc(ii-1)+maxIndc),singleCellTrace(peakLoc(ii-1)+maxIndc),'r.');
%                 plot(singleCellTraceTime(peakLoc(ii)+minIndc),singleCellTrace(peakLoc(ii)+minIndc),'r.');
                          
                temp_OffSlope(ii) = minOffSlope;
                temp_OffSlopeTime(ii) = singleCellTraceTime(peakLoc(ii)+minIndc);
                temp_OnSlope(ii) = maxOnSlope;
                temp_OnSlopeTime(ii) = singleCellTraceTime(peakLoc(ii-1)+maxIndc);
                temp_PulseWidth(ii) = singleCellTraceTime(peakLoc(ii)+minIndc) - singleCellTraceTime(peakLoc(ii-1)+maxIndc);
            end
            
        end
        
        if ii ~= 1 && ii+1 <= numPeaks

            % a middle peak 
            onSlope = diff(singleCellTrace(peakLoc(ii-1):peakLoc(ii)));
            [maxOnSlope, maxIndc]=max(onSlope);

            offSlope = diff(singleCellTrace(peakLoc(ii):peakLoc(ii+1)));
            [minOffSlope, minIndc]=min(offSlope);



            
            % if there is no OffSlope or PulseWidth associated with peak,
            % then put NaN. Disregard peak.
            if isempty(minOffSlope)|isempty(maxOnSlope)
                temp_OffSlope(ii) = NaN;
                temp_OffSlopeTime(ii) = NaN;
                temp_OnSlope(ii) = NaN;
                temp_OnSlopeTime(ii) = NaN;
                temp_PulseWidth(ii) = NaN;
            else
%                 % test
%                 figure(2); plot(singleCellTraceTime, singleCellTrace); hold on; 
%                 plot(singleCellTraceTime(peakLoc(ii-1)+maxIndc),singleCellTrace(peakLoc(ii-1)+maxIndc),'r.');
%                 plot(singleCellTraceTime(peakLoc(ii)+minIndc),singleCellTrace(peakLoc(ii)+minIndc),'r.');
            
                temp_OffSlope(ii) = minOffSlope;
                temp_OffSlopeTime(ii) = singleCellTraceTime(peakLoc(ii)+minIndc);
                temp_OnSlope(ii) = maxOnSlope;
                temp_OnSlopeTime(ii) = singleCellTraceTime(peakLoc(ii-1)+maxIndc);
                temp_PulseWidth(ii) = singleCellTraceTime(peakLoc(ii)+minIndc) - singleCellTraceTime(peakLoc(ii-1)+maxIndc);
            end
            
        end
        
    end
    
    Features.OnSlope = temp_OnSlope;
    Features.OnSlopeTime = temp_OnSlopeTime;
    Features.OffSlope = temp_OffSlope;
    Features.OffSlopeTime = temp_OffSlopeTime;
    Features.PulseWidth = temp_PulseWidth;
    
    clear temp_OnSlope; clear temp_OnSlopeTime; clear temp_OffSlope; clear temp_OffSlopeTime;
    clear temp_PulseWidth;
    
end
end





%    'height_distr', 'width_distr', 'time_to_peak_distr', 'autocorr_distr')