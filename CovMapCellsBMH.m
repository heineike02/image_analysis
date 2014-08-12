function timecoursedata =CovMapCellsBMH(imdir,images, channel, circ, imbg , siz , storeim, rad, std_thresh, time_vals)
%function which finds cells in a  given image
%di = folder containing .tif images (string), 
%images = list of image filenames to use in the form of a cell with N string entries.
%If there are M different channels, the cell is an NxM. 
%channel = string of channel being used.  Currently only one channel at a
%time
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
%path(path,'C:\Users\Ben\Google Drive\HES_Lab_Share\Scripts\JSO_Image_Analysis');

%N = length(images.(channels{1}));
N = length({images.(channel)});
if N ~= length(time_vals)
    'Error - times not the same size as number of images'
    return
end


%Nchan = length(channels);
% kk = 1; 
% for jj = 1:Nchan
%     if strfind(channels{jj}, 'BF')>0
%         'Bright Field images will be used'
%         bf_flag = 1;
%     else 
%         bf_flag = 0;
%         channels_epi{kk} = channels{jj};
%         kk = kk + 1;
%     end
% end


timecoursedata = struct('name','','celldata',[],'time',0);
%timecoursedata(1:N,1:(Nchan-bf_flag)) = timecoursedata(1);
timecoursedata(1:N,1) = timecoursedata(1);

for jj = 1:N
    imnames = {images.(channel)};
    imname = imnames(jj);
    strcat(int2str(jj), {' of '}, int2str(N), {' Channel = '}, channel, {' '}, imname  )
    im = imread(char(strcat(imdir,imname)));
%     if bf_flag == 1
%         bf_mask = imread(char(strcat(imdir,images.BF{jj})));
%    else
    bf_mask = [];
%    end
    celldata = FindCellsBMH(im, circ, bf_mask , imbg , siz, storeim, rad, std_thresh);
    %change this if I ever add another channel
    kk = 1
    timecoursedata(jj,kk).name = imname;
    timecoursedata(jj,kk).celldata = celldata;
    timecoursedata(jj,kk).time = time_vals(jj);
end


end

