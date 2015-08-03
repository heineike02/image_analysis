function timecoursedata = CovMapCells(imdir, images, time_vals, channels_to_image, cell_find_params)
%function which finds cells in a set of images
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

imbg = cell_find_params.imbg;

N = length(images.(channels_to_image{1}));

if N ~= length(time_vals)
    'Warning: times not the same size as number of images.  This will occur if you removed bad images from the source directory and the metadata no longer matches the number of files'
    %return
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

    for ch = 1:length(channels_to_image);
        channel = channels_to_image{ch};
        imnames = images.(channel);
        imname = imnames(jj);
        cell_find_params.imbg = imbg.(channel);
        strcat(int2str(jj), {' of '}, int2str(N), {' Channel = '}, channel, {' '}, imname  )
        im = imread(char(strcat(imdir,imname)));
        celldata.(channel) = FindCells(im, cell_find_params);

    end

    %combine celldata from each channel into one
    L = cell_find_params.rad;
    %celldata = match_two_channels(celldataCH.(channels{1}),channels{1},celldataCH.(channels{2}),channels{2},L);
    %My algorithm breaks on 20140703 A11_site1.  
    if length(channels_to_image)==2
        celldata = match_two_channels_Hung(celldata.(channels_to_image{1}),channels_to_image{1},celldata.(channels_to_image{2}),channels_to_image{2},L);
        %Make this function handle a situation where no cells are detected.
    elseif length(channels_to_image)==1
        %makes celldata have a channel field if there is only one channel
        %to remain consistant with multiple channel image data structure.
        channel = channels_to_image{1};
        celldata_new = struct();
        for mm = 1:length(celldata.(channel));
            celldata_fields = fields(celldata.(channel));
                for kk = 1:length(celldata_fields)
                    celldata_new(mm).(celldata_fields{kk}).(channel) = celldata.(channel)(mm).(celldata_fields{kk});
                end
            celldata_new(mm).Cxloc.Combined = celldata_new(mm).Cxloc.(channel); 
            celldata_new(mm).Cyloc.Combined = celldata_new(mm).Cyloc.(channel);
        end 
        celldata = celldata_new;
    else
        error('NchanError','more than three channels not yet supported')
    end
    

%     else
%        imnames = images.(channel);
%        imname = imnames(jj);
%        strcat(int2str(jj), {' of '}, int2str(N), {' Channel = '}, channel, {' '}, imname  )
%        im = imread(char(strcat(imdir,imname)));
%        cell_find_params.imbg = imbg.(channel);
%        celldata = FindCells(im, cell_find_params);
% 
%     end

    timecoursedata(jj).name = imname;
    timecoursedata(jj).celldata = celldata;
    clear celldata
    timecoursedata(jj).time = time_vals(jj);
end

end


