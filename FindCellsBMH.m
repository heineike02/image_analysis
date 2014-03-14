function celldata = FindCellsBMH(im, circ, bf_mask, imbg , siz, storeim, rad, std_thresh)
%Note: imbg should be median filtered and converted to double to avoid
%artifacts.

%number of pixels that smoothing uses
coarse_smooth = 25;
local_smooth = 3; 

max_cells = 500;
%number greater than the maximum number of cells in an image

%Number of iterations to use in the deconvlucy algorithm.
iterations = 5;

%Parameters for finding maxima
close_max = 1;
far_max = 6;

%Number of pixels to use in calculating nuclear enrichment - 
ne_pixels = 5;
%nuclear enrichment is defined as the ratio of the intensity of the
%top five pixels within the rad distance of the center of the image (where the 
%image is brightest) 
%and the median intensity of the pixels within the rad distance.  


%currently not capable of analyzing two channels.  
%If I do implement two channels, I want to have the images as 3d arrays -
%2d for the images and one dimension for channels.  
%also for two channels if the location of a cell is off, how would you
%choose between the two?

%adds JSO image processing files to path
%Windows
%path(path,'C:\Users\Ben\Google Drive\HES_Lab_Share\Scripts\JSO_Image_Analysis');

%function which finds cells in a given image
%circ=the image of the platonic ideal of a cell (2d matrxi of doubles), 
%bf_mask = file name of bright field mask
%storim=flag(0=only stores summary statistics, 1=store image of each cell),
%siz(1) is the number of square pixels used to store the image of the largest cell in the
%image.  siz(2) is the number of square pixels used to calculate nuclear enrichment. The default value for both is 17. 
%rad = cell's radius for nuclear enrichment.  Current typical cell radius is 6

%Load the default background image: 
%          imbg=double(imread('RFP_BG.tif'));

%Load the default circ
%          load ../circ.mat


%metadata = file name of metadata file


%output: 
%N - number of cells
%
%celldata array of structures for each cell with fields:
%image: double image
%nf : nuclear enrichment
%nmi : median brightness of each cell
%Cxloc
%Cyloc




%Want to have option of using masking BF image
%should probably take a background image to normalize for variation each
%time

%Sets Defaults    

%Default is to store images, to use the default background image 'RFP_BG.tif', and to not use a bright field mask
if nargin<6
    storeim=1   
    if nargin<5
        siz = [17,17]
        if nargin<4
          %if no background image is entered, imbg is equal to 1
          imbg = 1
          %kr=0; qr=0;
          if nargin<3
             bf_mask = 0
          end
        end
    end

end
%     

%Block size is a certain number of pixels away from the center in either
%direction.  Called a radius even though it is a square block
block_radius = floor(siz(1)/2);

%build arrays to store output


%%%Not sure if this should be a double
%im_init = uint16(0);
%im_init(siz,siz) = 0;
im_init = zeros(siz(1),siz(1));

if storeim == 1 
    celldata = struct('image',im_init,'nf',0.0,'nmi',0.0,'Cxloc',0.0,'Cyloc',0.0);
else
    celldata = struct('nf',0.0,'nmi',0.0,'Cxloc',0.0,'Cyloc',0.0);
end

celldata(1:max_cells) = celldata(1);

%smothkern=[0.3,0.3,0.3;0.3,1,0.3;0.3,0.3,0.3];
%Smoothing kernel is 3x3 100% in the middle, and 10% of surrounding cells.
%Doesn't look like the smoothing kernel is ever used...
smothkern=[0.1,0.1,0.1;0.1,1,0.1;0.1,0.1,0.1];
     
%load first image
%im=double(imread([di,'\',char(files(in(i)).name)]));
im = double(im);
%normalize to background (be careful here, not always a good idea)
im=im./imbg;

%deconvolute with kernal of ideal cell, to find cells as peaks

ima=deconvlucy(im,circ,iterations);
        
%smooth image heavily
imcs=(ima./conv2(ima,ones(coarse_smooth),'same'));
%smooth locally
%%%originally this one didn't have 'same' - not sure why
imc=conv2(imcs,ones(local_smooth),'same');
%find possible locations of cells with generous threshold
%threshold is median plus twice the mean absolute deviation. 
thr=median(imc(:))+2*mad(imc(:));

%Not sure what thr_map is for
thr_map = imc>thr;
%assumes square images

imsize = size(im);

%strange condition - don't quite understand it - why does it just sum along
%one dimension??
%resets threshold to the 90th quantile if there are more than 10% pixels greater than the threshold in one column of the image.  
%if sum(thr_map)/(imsize(1)*imsize(2)) > 0.10 ; thr=quantile(ima(:),0.90); 'new threshold'; end
        
%find possible cell pixels
[px,py]=find(imc>thr);

%For troubleshooting
%plot(py,px,'rx')

%find local maxima to define cells
%arguments: 
%im - image to analyze
%px - x coordinates of pixels above threshold
%py - y coordinates of pixels above threshold
%dis - distance (in pixels) around which max will be defined

%not sure if the center of the cell should be found with imc instead of
%ima.

[maxx,maxy]=FindMaximaBMH(ima,px,py,close_max);
[maxx,maxy]=FindMaximaBMH_2(ima,maxx,maxy,far_max);        




% if yinc
%     cffT=zeros(2000,2*(siz(2)^2));
% else
%     cffT=zeros(2000,(siz(2)^2));
% end

%Sets up vectors to store cell coordinates
CxlocT = zeros(max_cells,1);
CylocT = zeros(max_cells,1);

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
           imb=GetBlock(im,block_radius,xloc,yloc);
           %only keeps image if the standard deviation / mean of the pixels
           %is greater than std_thresh
           %std_thresh = 0.2;
           if std(imb(:))/mean(imb(:))>std_thresh
            %adjusts the location of the center using the single cell image
            %Not sure what the mvt was in the original code - just set to
            %one in the function.  I took that out. 
            %mvt = 1;
            xy = lcrosscorr(imb);
            yloc=yloc+(xy(2)-block_radius);
            xloc=xloc+(xy(1)-block_radius);
            %check to see if cell with new center is too close to the edge
            if (xloc + block_radius) <= imsize(1) & (xloc-block_radius) > 0 & (yloc + block_radius) <= imsize(2) & (yloc-block_radius)>0
                %check to see proposed new center is already on the list
                if sum((xvec == xloc) & (yvec==yloc)) == 0
                    qr=qr+1;
                    if qr >= max_cells
                        ['more than max number of cells ', num2str(max_cells)];
                    end
                    %mv(qr)=mvt;  %seems to be 1 no matter what don't see why
                    %we save this.  
                    im_out = GetBlock(im,block_radius,xloc,yloc);
                    if storeim == 1
                        celldata(qr).image = im_out;
                    end
                    celldata(qr).Cxloc = xloc;
                    celldata(qr).Cyloc = yloc;
                    xvec = [xvec,xloc];
                    yvec = [yvec,yloc];
                    [nf,nmi]= NEnrichBMH(im_out,circ, rad, ne_pixels);
                    celldata(qr).nf = nf;
                    celldata(qr).nmi = nmi;
                    %if yinc; [nf,nmi]=NEnrich(reshape(cffT(j,(siz(2)^2)+1:2*(siz(2)^2)),siz(2),siz(2)),circ); cff(kr+j,3:4)=[nf,nmi]; end
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


function [acqd,strain,numw,numc]=readAcqDat(fi,sfi)
%open acqdat.txt
fid=fopen(fi);
i=0;
while ~feof(fid)
   i=i+1;
   l=fgetl(fid);
   
   c=(regexpi(l,'\t'));
   el=length(l);
   
   t(i)=str2num(l(1:c(1)-1));
   
   k=1;
   wellnum(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   LocX(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   LocY(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   LocZ(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   PFoff(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   ch(i)=cellstr(l(c(k)+1:c(k+1)-1)); k=k+1;
   exp(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   
   if length(c)>9
       ledsetval(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
       ledoutval(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
       fname(i)=cellstr(l(c(k)+1:el(1))); 
   else
       fname(i)=cellstr(l(c(k)+1:el(1)));  
   end
end

fclose(fid);

numw=length(unique(wellnum));
numc=length(unique(ch));    
acqd.wellnum=wellnum;
acqd.ch=ch;
acqd.time=t;
acqd.exp=exp;
if length(c)>9
    acqd.LEDset=ledsetval;
    acqd.LEDout=ledoutval;
end
acqd.fname=fname;

i=0;
if ~isempty(sfi)
    %open strains.txt
    fid=fopen(sfi);
    while ~feof(fid)
        i=i+1;
        strain(i)=cellstr(fgetl(fid));
    end
    fclose(fid)
else
    strain(1:(numw))=cellstr(numw);
end

end


function [cond]=readconDat(fi)
%read conditions.txt return array

i=0;
    fid=fopen(fi);
    while ~feof(fid)
        i=i+1;
        cond(i)=cellstr(fgetl(fid));
    end
    fclose(fid)
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
