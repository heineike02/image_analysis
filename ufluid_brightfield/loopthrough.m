function TIMECOURSEDATA = loopthrough(IMDIR,IMAGES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION - generates a storage structure for each image frame and for    %
%each cell                                                                %
%                                                                         %
%INPUTS:                                                                  %
%           IMDIR -  folder containing .tif images (string)               %
%           IMAGES - image filenames list in the form of a cell with N    %
%                    string entries.                                      %
%           CHANNELS - cell with strings containing image channels.       %
%                      Ex: {'RFP','YFP'}                                  %
%                                                                         %
%OUTPUTs:                                                                 %
%           TIMECOURSEDATA - structure with fields                        %
%                                   name    = filename                    %
%                                   stats   = stats structure from        %
%                                             bfanalyze                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%N = number of images to loop through
N = length(IMAGES);

%initialize structure to store all information
TIMECOURSEDATA = struct('name','','stats',[]);

%loop through each frame and find cells + relevant metrics
for jj = 1:N 
        %display progress
        [int2str(jj), ' of ', int2str(N), ' ', IMAGES(jj).FLUOR]
        
        %generate input images
        BFim = imread([IMDIR,IMAGES(jj).BF]);
        FLUORim = imread([IMDIR, IMAGES(jj).FLUOR]);
        
        %determine whether to translate in XY direction in fluorescence im
        if jj == 1
            shift = 1;
        else
            shift = 0;
        end
        
        %cell segmentation and metrics
        STATS = bfanalyze(BFim, FLUORim, shift);
        
        %populate storage structure
        TIMECOURSEDATA(jj).name = IMAGES(jj).FLUOR;
        TIMECOURSEDATA(jj).stats = STATS;
end

