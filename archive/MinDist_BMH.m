function timecoursedata_out = MinDist_BMH(di, timecoursedata , max_shift_rad, min_dist_thresh, storeim) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from MinDist_SYC to accommodate output from CovMapCellsBMH
%
% 2. Adjust xy locations on second image accordingly
% 3. Perform MinDist Calculation

%di = folder containing .tif images (string), 
%timecoursedata - structure with fields
%name = filename
%celldata = celldata structure from FindCellsBMH - length is number of
%cells identified in that frame.
%time = time at which image was taken. 
%
%circ=the image of the platonic ideal of a cell (2d matrxi of doubles),
%max_shift_rad = max distance (in pixels) we expect the image to shift
%between frames. 
%timecoursedata_out = structure with same fields, but indices of cells
%correspond to indices from the first frame.  

timecoursedata_out(1) = timecoursedata(1);
imdata_old = timecoursedata(1).celldata;
im_old = imread([di,timecoursedata(1).name]);

%for each image in the timecourse after the initial image
for jj = 2:length(timecoursedata);
    [int2str(jj), ' of ', int2str(length(timecoursedata)),' Cell Tracking']
    imdata_new = timecoursedata(jj).celldata;
    im_new = imread([di,timecoursedata(jj).name]);
    
    %initialize output value
    clear('imdata_out')
    if storeim == 1 
        imdata_out = struct('image',zeros(size(im_old)),'nf',0.0,'nmi',0.0,'Cxloc',0.0,'Cyloc',0.0);
    else
        imdata_out = struct('nf',0.0,'nmi',0.0,'Cxloc',0.0,'Cyloc',0.0);
    end
    imdata_out(1:length(imdata_old)) = imdata_out(1);
    
    x_old = [imdata_old.Cxloc];
    y_old = [imdata_old.Cyloc];
    x_new = [imdata_new.Cxloc];
    y_new = [imdata_new.Cyloc];
    
    % Calculate shift between two images 
    xy_shift = ImageShift_BMH(im_old, im_new, max_shift_rad);
    
    % Adjust xy locations on second image
    x_new = x_new + xy_shift(1);
    y_new = y_new + xy_shift(2);
    
    % for each cell in the old image   
    for kk = 1:length(imdata_old);
        xx = x_old(kk);
        yy = y_old(kk);
        dist = zeros(length(imdata_new),1);
        %if cell in old image is NaN, assign NaN set to output image
        if imdata_old(kk).Cxloc == NaN;
            imdata_out(kk).Cxloc = NaN;
            imdata_out(kk).Cyloc = NaN;
        else
            %Find minimum distance between centers
            for mm = 1:length(imdata_new);
                xx2 = x_new(mm);
                yy2 = y_new(mm);
                dist(mm) = sqrt((xx-xx2)^2+(yy-yy2)^2);
            end    
            [min_dist,min_ind] = min(dist);
            %If minimum distance is below a threshold
            if min_dist < min_dist_thresh;
                %Assign cell in output image data to current cell from old image
                imdata_out(kk) = imdata_new(min_ind);
                %remove assigned cell from new image and x/y locations
                imdata_new(min_ind) = [];
                x_new(min_ind) = [];
                y_new(min_ind) = [];
            else
                %else assign an NaN set to that place in output image
                %data
                imdata_out(kk).Cxloc = NaN;
                imdata_out(kk).Cyloc = NaN;
            end
            
        end
    end
    
    timecoursedata_out(jj) = timecoursedata(jj);
    timecoursedata_out(jj).celldata = imdata_out;
    imdata_old = imdata_out;
    im_old = im_new;
end
      
    


%metrics for determining cell alignment
%1) distance from previous frame
%2) number of neighbors as compared to previous frame
%implicit assumption is 1st frame is reference frame for rest of frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %find all nuclear intensities of every single cell in every frame in the
% %same order as cellx_new and celly_new
% [cellV]= NucIntensity(cff, tim, Cyloc, Cxloc, circ);
% 
% %initialize cellx and celly and cellV
% cellx = nan(max(tim), 1000);
% celly = nan(max(tim), 1000);
% cellV_init = nan(max(tim), 1000);
% %%
% %establish 1st frame
% %indices for first frame to access the x and y locations in the first frame
% indx_1 = find(tim==1);
% indy_1 = find(tim==1);
% %finding the xy locations of the first frame
% xloc_1 = Cxloc(indx_1);
% yloc_1 = Cyloc(indy_1);
% %fill the first row of the output xy locations
% %note that these initial seed cells in frame one will be the standard
% %numbering that all other frames use
% cellx(1,1:length(xloc_1)) = xloc_1;
% celly(1,1:length(yloc_1)) = yloc_1;
% cellV_init(1,1:length(cellV(1).cellnum)) = cellV(1).cellnum;
% %note the index of this first frame
% 
% %new version of cellx and celly
% cellx_new = cellx(:,find(~isnan(cellx(1,:))));
% celly_new = celly(:,find(~isnan(celly(1,:)))); 
% cellV_new = cellV_init(:, find(~isnan(cellV_init(1,:))));
% 
% %call neighborhood algorithm to find neighborhood cells for all frames and
% %all cells
% [cellstruct] = Neighborhood_SYC(Cxloc, Cyloc, tim);

% %readout for all the cells
% im = imread(strcat(di, 'RFP_p1_t1.tiff'));
% figure(1);
% imagesc(im);
% hold on;
% H= plot(celly(1,1:length(yloc_1)), cellx(1,1:length(xloc_1)), 'r.');
% set(H, 'MarkerSize', 8);

% %readout just one cell 
% %with index in celly and cellx of 1
% figure(1);
% im = imread(strcat(di, 'RFP_p1_t1.tiff'));
% imagesc(im);
% hold on;
% H= plot(celly_new(1,10), cellx_new(1,10), 'r.');
% set(H, 'MarkerSize', 8);

%%
%loop through all cells in 1st frame

% for j = 1:length(timecoursedata(1).)
% %trace the same cell as the cell in frame 1
% %loop through all cells in each frame
% for i = 2:length(celly_new(:,j))
%     x = Cxloc(find(tim == i));
%     y = Cyloc(find(tim == i));
%     cellV_frame = cellV(i).cellnum;
%     
%     %distance criterion
%     dx = x - cellx_new(1,j);
%     dy = y - celly_new(1,j);
%     ddist = sqrt(dx.^2+ dy.^2);
%     
%     [dist_val dist_ind] = sort(ddist);
%     
%     %take top5 minimum distance candidates
%     for k = 1:5
%         bin_stor = zeros(1,5); %initiate a binary storage vector
%         
%         %compare first frame neighbors to current frame neighbors
%         neigh1 = cellstruct(1).cells(j).len; %first frame number of neighbors
%         neigh2 = cellstruct(i).cells(dist_ind(k)).len; %current frame number of neighbors
%         if neigh1 == neigh2
%             %fprintf('same number of neighbors k=%d', k)
%             %check that relative shifts are the same
%             %x location check
%             currx = x(dist_ind(k)); %current frame xlocation
%             firx = cellx_new(1,j); %first frame xlocation
%             %y location check
%             curry = y(dist_ind(k)); %current frame ylocation
%             firy = celly_new(1,j); %first frame ylocation
%             
%             if neigh1 ~= 0 
%                 neigh_currx = cellstruct(i).cells(dist_ind(k)).x; %neighbors current xlocation
%                 neigh_firx = cellstruct(1).cells(j).x; %neighbors first xlocation
%             
%                 relx = abs(currx - firx); %difference 
%                 neigh_relx = abs(neigh_currx - neigh_firx); %vector of some length containing differences 
%                 mean_neigh_relx = mean(neigh_relx);
%            
%                 neigh_curry = cellstruct(i).cells(dist_ind(k)).y; %neighbors current ylocation
%                 neigh_firy = cellstruct(1).cells(j).y; %neighbors first ylocation
%             
%                 rely = abs(curry - firy); %differences
%                 neigh_rely = abs(neigh_curry - neigh_firy); %vector of some length containing differences
%                 mean_neigh_rely = mean(neigh_rely);
%             
%                 %relative shifts from current frame to 1st frame the difference
%                 %between these numbers should be small
%                 reldist = sqrt(relx.^2+rely.^2);
%                 neigh_reldist = sqrt(mean_neigh_relx.^2+ mean_neigh_rely.^2);
%             
%                 %calculate relative shift - smaller the better
%                 relshift = (reldist - neigh_reldist).^2;
% 
%                 bin_stor(k) = relshift;
%             else
%                 bin_stor(k) = NaN; %this can be improved - by using a fixed cell as a reference point and taking the relshift
%             end
%             
%         else
%             bin_stor(k) = NaN; %just a large number 
%         end 
%         [minval, minind] = min(bin_stor);
%     end
%     
%     %neighborhood option
%     %assign the cell in the current frame to a cell in the 1st frame
%     cellx_new(i,j) = x(dist_ind(minind));
%     celly_new(i,j) = y(dist_ind(minind));
%     
%     cellV_new(i,j) = cellV_frame(dist_ind(minind));
%     %fprintf('i=%d\n', i)
%     
%     %distance but no neighborhood option 
%     %cellx_new(i,j) = x(dist_ind(1));
%     %celly_new(i,j) = y(dist_ind(1));
%     
%     [rangex, rangey, cellx_newer, celly_newer, cellV_newer] = straightcells(cellx_new, celly_new, cellV_new);
% end
% 
% end

end
 