function celldata = FindCells(im, cell_find_params)
%Note: imbg should be median filtered and converted to double to avoid
%artifacts.

circ = cell_find_params.circ;
imbg = cell_find_params.imbg;
siz = cell_find_params.siz;
storeim = cell_find_params.storeim;
rad = cell_find_params.rad;
std_thresh = cell_find_params.std_thresh;
thresh = cell_find_params.thresh;

%number of pixels that smoothing uses
coarse_smooth = cell_find_params.coarse_smooth; 
local_smooth = cell_find_params.local_smooth; 

maxcells = cell_find_params.maxcells;
%number greater than the maximum number of cells in an image

%Number of iterations to use in the deconvlucy algorithm.
deconvlucy_iterations = cell_find_params.deconvlucy_iterations;

%Parameters for finding maxima
close_max = cell_find_params.close_max;
far_max = cell_find_params.far_max;
edge_margin = cell_find_params.edge_margin; 

%Number of pixels to use in calculating nuclear enrichment - 
ne_pixels = cell_find_params.ne_pixels;
%nuclear enrichment is defined as the ratio of the intensity of the
%top five pixels within the rad distance of the center of the image (where the 
%image is brightest) 
%and the median intensity of the pixels within the rad distance.  


%function which finds cells in a given image
%circ=the image of the platonic ideal of a cell (2d matrxi of doubles), 
%bf_mask = file name of bright field mask
%storim=flag(0=only stores summary statistics, 1=store image of each cell),
%siz(1) is the number of square pixels used to store the image of the largest cell in the
%image.  siz(2) is the number of square pixels used to calculate nuclear enrichment. The default value for both is 17. 
%rad = cell's radius for nuclear enrichment.  Current typical cell radius is 6


%Block size is a certain number of pixels away from the center in either
%direction.  Called a radius even though it is a square block
block_radius = floor(siz(1)/2);

im_init = zeros(siz(1),siz(1));

if storeim == 1 
    celldata = struct('image',im_init,'nf',0.0,'nmi',0.0,'Cxloc',0.0,'Cyloc',0.0);
else
    celldata = struct('nf',0.0,'nmi',0.0,'Cxloc',0.0,'Cyloc',0.0);
end

celldata(1:maxcells) = celldata(1);
     
%load first image
%im=double(imread([di,'\',char(files(in(i)).name)]));
im = double(im);
%normalize to background (be careful here, not always a good idea)
im=im./imbg;

%deconvolute with kernal of ideal cell, to find cells as peaks DEBUG POINT
%imagesc(ima)
ima=deconvlucy(im,circ,deconvlucy_iterations);
        
%smooth image heavily
imcs=(ima./conv2(ima,ones(coarse_smooth),'same'));
%smooth locally
%%%originally this one didn't have 'same' - not sure why DEBUG POINT
%imagesc(imc)
imc=conv2(imcs,ones(local_smooth),'same');
%find possible locations of cells with generous threshold
%threshold is median plus thresh* the mean absolute deviation. 
thr=median(imc(:))+thresh*mad(imc(:));

imsize = size(im);
       
%find possible cell pixels
[px,py]=find(imc>thr);

%DEBUG POINT For troubleshooting all pixels above threshold
%plot(py,px,'rx')

%find local maxima to define cells
%arguments: 
%im - image to analyze
%px - x coordinates of pixels above threshold
%py - y coordinates of pixels above threshold
%dis - distance (in pixels) around which max will be defined

%not sure if the center of the cell should be found with imc instead of
%ima.


%run once with close max not removing duplicates.
[maxx,maxy]=FindMaxima(ima,px,py,close_max,0,edge_margin);       
%run again with far max removing duplicates. 
[maxx,maxy]=FindMaxima(ima,maxx,maxy,far_max,1,edge_margin);
% DEBUG POINT - check if duplicates are removed plot(maxy,maxx,'rx')


%Sets up vectors to store cell coordinates
CxlocT = zeros(maxcells,1);
CylocT = zeros(maxcells,1);

%go through putative cells, pull pixel region, analyze them

%index for new vector of cells - will only keep cells if certain conditions
%are met
qr=0;
xvec = [];
yvec = [];
for j=1:length(maxx);
        yloc=maxy(j);
        xloc=maxx(j);
        %throws out cells that are too close to the edge
        if (xloc + block_radius) <= imsize(1) & (xloc-block_radius) > 0 & (yloc + block_radius) <= imsize(2) & (yloc-block_radius)>0 
           % DEBUG POINT - check to see if vaguely looks like cell (may be
            % off center, should be bright pixel in center) imagesc(imb)
            imb=GetBlock(im,block_radius,xloc,yloc);
           %only keeps image if the standard deviation / mean of the pixels
           %is greater than std_thresh
           % DEBUG POINT - check to see if cell being thrown out, is number
           % on left of > greater than std_thresh? NEED TO FIX - SHOULD BE
           % EVALUATING AFTER CENTERING IMAGE
           if std(imb(:))/mean(imb(:))>std_thresh
             
            %adjusts the location of the center using the single cell image
            xy = lcrosscorr(imb);
            yloc=yloc+(xy(2)-block_radius);
            xloc=xloc+(xy(1)-block_radius);
            %check to see if cell with new center is too close to the edge
            if (xloc + block_radius) <= imsize(1) & (xloc-block_radius) > 0 & (yloc + block_radius) <= imsize(2) & (yloc-block_radius)>0
                %check to see proposed new center is already on the list
                if sum((xvec == xloc) & (yvec==yloc)) == 0
                    qr=qr+1;
                    if qr >= maxcells
                        ['more than max number of cells ', num2str(maxcells)];
                    end
                    im_out = GetBlock(im,block_radius,xloc,yloc);
                    % DEBUG POINT Final picture of cell processed - check this to see
                    % if it looks like a cell. imagesc(im_out)
                    if storeim == 1
                        celldata(qr).image = im_out;
                    end
                    celldata(qr).Cxloc = xloc;
                    celldata(qr).Cyloc = yloc;
                    xvec = [xvec,xloc];
                    yvec = [yvec,yloc];
                    [nf,nmi]= NEnrich(im_out,circ, rad, ne_pixels);
                    celldata(qr).nf = nf;
                    celldata(qr).nmi = nmi;                    
                end
                
            end
           end
        end
end   

%Shrink down celldata to fit number of cells
celldata = celldata(1:qr);
end


    
function im_block = GetBlock(im,block_radius,locX,locY)
    %load a portion of an image
    %not set up for channel data.  To do channel data would want to have im
    %contain all channels of the image in different dimensions. 
    
    %Not sure why we are passing around images as long vectors. removed
    %that.  
    
    im_block = im(locX-block_radius:locX+block_radius,locY-block_radius:locY+block_radius);
    
%     
%     if ~isempty(im2)
%         rim=im2(locX-hd:locX+hd,locY-hd:locY+hd);
%         bl=[tim(:)',rim(:)'];
%     else
%         bl=tim(:)';
%     end

end


function xy = lcrosscorr(im)
    %gives the max of a convolution with ones of the image. 
    %convolution with ones seems to smooth the image.
    %Problem: if there is two cells in the image, this only finds the
    %center of one. could we make it find multiple local maxima and only
    %pick the max that is closest to the center? 
    imc=conv2(im,ones(5),'same');
    [a1,row_vec]=max(imc);
    [a2,col] = max(a1);
    xy = [row_vec(col),col];
    
end
