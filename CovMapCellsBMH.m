function timecoursedata =CovMapCellsBMH20130716(di,images, channels, circ, images_bf, imbg , siz , storeim, rad, std_thresh, times)

%function which finds cells in a  given image
%di = folder containing .tif images (string), 
%images = list of image filenames to use in the form of a cell with N string entries.
%If there are M different channels, the cell is an NxM. 
%channels = cell with strings containing image channels.  Ex: {'RFP','YFP'} 
%circ=the image of the platonic ideal of a cell (2d matrxi of doubles), 
%bf_images = bright field images should be the same length as images.
%imbg = background images - should be one for each channel??
%siz= a 2by1 vector of doubles, 
%storim = flag(0=only stores summary statistics, 1=store image of each cell), 
%times = vector of times - should be the same length as images (one time
%for each image).
%
%output: 
%timecoursedata - structure with fields
%name = filename
%celldata = celldata structure from FindCellsBMH - length is number of
%cells identified in that frame.
%time = time at which image was taken.  
%
%adds JSO image processing files to path
%Windows
path(path,'C:\Users\Ben\Google Drive\HES_Lab_Share\Scripts\JSO_Image_Analysis');

N = length(images);
if N ~= length(times)
    'Error - times not the same size as number of images'
    return
end
    
timecoursedata = struct('name','','celldata',[],'time',0);
timecoursedata(1:N) = timecoursedata(1);

for jj = 1:N
    [int2str(jj), ' of ', int2str(N)]
    im = imread([di,images{jj}]);
    bf_mask = imread([di,images_bf{jj}]);
    celldata = FindCellsBMH(im, circ, bf_mask , imbg , siz, storeim, rad, std_thresh);
    timecoursedata(jj).name = images{jj};
    timecoursedata(jj).celldata = celldata;
    timecoursedata(jj).time = times(jj);
end


end

