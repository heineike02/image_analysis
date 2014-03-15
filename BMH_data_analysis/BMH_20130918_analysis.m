% Incorporate Cell Tracking
% Try with K. Lactis
% Visualize individual and mean nuclear localization for glucose shift
% before and after Glucose deprivation in S.C. and K. Lactis
% Repeat for BPAC (in another file)

%load image - first image is S. Cerevisiae pre-induction
%Windows: 
ipdir = 'C:\Users\Ben\Google Drive\HES_Lab_Share\Scripts\JSO_Image_Analysis\BMH_Image_Analysis\'
imdirSP.KL = 'C:\Users\Ben\Documents\Data\PKA_Project\20130918\Blue_light_Strain_K_1\'
imdirSP.SC = 'C:\Users\Ben\Documents\Data\PKA_Project\20130918\Blue_light_Strain_SC_BPAC\'

%Mac
%ex: imdir = '/Users/ucsf/Google Drive/UCSF/ElSamad_Lab/PKA/WetLab/20130716/Gluc_pre_ind/Pos1/'
storeim = 1;
imbg = 1;
%load image of KLactis cell - this was made in BMH_20130916_analysis_klac.m
load([ipdir, 'circKL.mat'])

%load image of S.Cerevisiae cell (made by Jacob)
load([ipdir(1:end-19),'circ.mat'])

%[celldata,N] = FindCellsBMH(im, circ, bf_mask, im_bg , siz, storeim)

%Size of image of individual cell to analyze
sizSP.KL = [15,15];
sizSP.SC = [17,17];
%radius of cell for nuclear enrichment calculation
radSP.KL = 4;
radSP.SC = 6;
%std thresh for calling a peak to find a cell
%0.16 seemed to work best for k.lactis on 16JUL images
std_threshSP.KL = 0.16;
std_threshSP.SC = 0.2;

max_shift_rad = 4;
min_dist_thresh = 3;

%There are four positions per condition.  

species_cell = {'KL'}
channels = {'RFP', 'YFP'}
%Eventually should track cells using both channels, but for now RFP is the
%best channel
cell_track_chan_ind = 1
Nchan = length(channels)
positions = 4
dt = .5


for sp = 1:length(species_cell)
    species = species_cell{sp};
    siz = sizSP.(species);
    rad = radSP.(species);
    std_thresh = std_threshSP.(species);
    imdir = imdirSP.(species);
    D = dir(imdir);
    %subtract 6 for the non-image files in the directory.
    Ntimes = (length(D)-6)/(Nchan*positions);
    times = dt*[0:Ntimes-1];

    %initializing plotting variables
    nf_vec_mean = zeros([Ntimes,positions]);
    
    for nn = 1:positions
        pos = num2str(nn);
        for mm = 1:Nchan
           images.(channels{mm}) = cell([Ntimes,1]);        
           for jj = 1:length(D)
               fname = D(jj).name;
               %check correct channel and position
               if strfind(fname,[channels{mm},'_p',pos]) > 0
                  %find time point and convert to index.  Channel ordering has to be in correct order (channel one t1, channel two t2 etc: 
                  t_ind = regexp(fname, 't(?<time>\d+)','names');
                  t_ind = str2double(t_ind.time);
                  kk = (t_ind + Nchan-mm)/Nchan;
                  images.(channels{mm}){kk} = D(jj).name;
               end
           end
        end
        
        %Get data for each position
        timecoursedata = CovMapCellsBMH(imdir ,images, channels, circ, imbg , siz , storeim, rad, std_thresh, times);
        %Just using RFP data to track cells
        timecoursedata_RFP = timecoursedata(:,cell_track_chan_ind);
        timecoursedata_out = MinDist_BMH(imdir, timecoursedata_RFP ,max_shift_rad, min_dist_thresh, storeim);
                
        %Plot data for each position in a subplot
        %count number of NaNs
  
        x_vec = zeros(length(timecoursedata_out(1).celldata),length(timecoursedata_out));
        y_vec = x_vec;
        nf_vec = x_vec;
        for jj = 1:length(timecoursedata_out);
            x_vec(:,jj) = [timecoursedata_out(jj).celldata.Cxloc];
            y_vec(:,jj) = [timecoursedata_out(jj).celldata.Cyloc];
            nf_vec(:,jj) = [timecoursedata_out(jj).celldata.nf];
        end
  
        figure(1)
        subplot(2,2,nn)
        plot(x_vec',y_vec')
        title('Position of cells')

        figure(2)
        subplot(2,2,nn)
        clf 
        hold on
        plot(times,nf_vec)
        xlabel('time (s)')
        ylabel('Nuclear Localization')
        title('Nuclear localization - Unfiltered')

        figure(3)
        subplot(2,2,nn)
        clf
        hold on
        nan_ind = ~isnan(x_vec(:,end));
        nf_vec_filt = nf_vec(nan_ind,:);
        plot(times,nf_vec_filt);
        title('Nuclear localization - filtered')
        xlabel('time (s)')
        ylabel('Nuclear Localization')

        figure(4)
        subplot(2,2,nn)
        clf
        hold on
        errorbar(times,mean(nf_vec_filt),std(nf_vec_filt))
        title('Mean Nuclear Localization')
        xlabel('time (s)')
        ylabel('Nuclear Localization')

        nf_vec_mean(:,nn) = mean(nf_vec_filt);

        
    end

      
end





%bf_mask = 'img_000000000_Default_000.tif';
%metadata = 'metadata.txt'

%S=load('/Users/ucsf/Google Drive/HES_Lab_Share/Scripts/JSO_Image_Analysis/circ.mat')
%circ=S.circ;
%circ=circ(1:13,2:14);

%}

