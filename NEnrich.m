function [nf,nmi]=NEnrich(im,circ, rad, ne_pixels)
%compute nuclear enrichment given an image and an ideal cell image
% Cell Radius: 
siz=size(im,1);
    %Why use a convolution (which reverses and shifts something) instead of
    %a cross-correlation which just shifts?
    a=conv2(double(im),double(circ),'same');
    
    [in,xy]=max(a(:));
    yloc=(mod(xy-1,siz(1))+1);
    xloc=floor(xy/siz(1))+1;
    
    for i=1:siz;
        matx(i,:)=1:siz;
        maty(:,i)=1:siz;
    end
    
    dis=sqrt((xloc-matx(:)).^2+(yloc-maty(:)).^2);
    
    %nuclear enrichment is defined as the ratio of the intensity of the
    %top five cells within the rad distance of the center of the image (where the 
    %image is brightest) 
    %and the median intensity of the pixels within the rad distance.  
    ne=sort(im(dis<rad),'descend');
    nmi = median(ne);
    nf=mean(ne(1:ne_pixels))./nmi;
end

    