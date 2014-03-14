function [cff,acqdat,Cxloc,Cyloc,tim]=CovMapCells(di,p,nframe,siz,circ,storim,mchan)

%function which find cells in image data
%di=folder containing .tif images (string), p=position number in image series (number),
%nfram=image numbers to analyze (array of doubles),siz= a 2by1 vector of
%doubles, circ=the image of the platonic ideal of a cell (2d matrxi of
%doubles), storim=flag(0=only stores summary statistics, 1=store image of
%each cell), mchan=first letter of the channel to be analyzed (eg.. RY is
%RFP than YFP)

if nargin<5 | isempty(circ);
    S=load('circ.mat')
    circ=S.circ;
    circ=circ(1:13,2:14);
end
if nargin<6
    storim=1;   
end
if nargin<7
    mchan='RY';       
end

%read the acqdat data and strains.txt txt file
[acqdat,strain,numw,numc]=readAcqDat([di,'\acqdat.txt'],[di,'\strains.txt']);


file = dir(strcat(di,'\*.tiff')); file2 = dir(strcat(di,'\*.tif'));
files=[file',file2'];

asiz=siz(1)^2;
k=0;

    for i=1:length(files);
        %process text of each file name
            n=regexpi(char(files(i).name),'_t(\d+).tif','tokens');
            nm(i)=str2num(char(n{1}));
            n=regexpi(char(files(i).name),'_p(\d+)_','tokens');
            np(i)=str2num(char(n{1}));
            n=regexpi(char(files(i).name),'(\w)FP_','tokens');
            ch(i)=(char(n{1}));
    end


    [x,in]=sort(nm);
    vin=in(find(np(in)==p & (ch(in)==mchan(2))));
    in=in(find(np(in)==p & (ch(in)==mchan(1))));
    hsiz=floor(siz/2);
    
    if isempty(vin); yinc=0;else yinc=1; end

    m=0;
    
    %load the backgound image
    imbg=double(imread('RFP_BG.tif'));
    kr=0; qr=0;
    
    %build arrays to store output
    if storim
        if yinc
            cff=zeros(50000,2*(siz(2)^2));
        else
            cff=zeros(50000,(siz(2)^2));
        end
    else
        if yinc
            cff=zeros(50000,4);
        else
            cff=zeros(50000,2);
        end        
    end
    
    if length(nframe)==1
        framestoload=1:nframe;
    else
        framestoload=nframe;
    end
    
    smothkern=[0.3,0.3,0.3;0.3,1,0.3;0.3,0.3,0.3];
     smothkern=[0.1,0.1,0.1;0.1,1,0.1;0.1,0.1,0.1];
     
     
     %loop through images
    for i=framestoload

        %load first image
        im=double(imread([di,'\',char(files(in(i)).name)]));
        
        %load second image
        if yinc
            im2=double(imread([di,'\',char(files(vin(i)).name)]));
            im2=im2./imbg;
        else
            im2=[];
        end
        
        %deconvolute with kernal of ideal cell, to find cells as peaks
        ima=deconvlucy((im)./imbg,circ,10);
        
        %smooth image heavily
        imc=(ima./conv2(ima,ones(25),'same'));
        %smooth locally
        imc=conv2(imc,ones(3));
        %find possible locations of cells with generous threshold
        thr=median(imc(:))+2*mad(imc(:));
         if sum(imc>thr)./(512^2)>0.10; thr=quantile(ima(:),0.90); end
        
         %find possible cell pixels
        [px,py]=find(imc>thr);
       
        %find local maxima to define cells
        [maxx,maxy]=FindMaxima2(ima,px,py,1);
        [maxx,maxy]=FindMaxima2(ima,maxx,maxy,6);        
        
        %normalize to background (be careful here, not always a good idea)
        im=im./imbg;
        
        qr=0;
        if yinc
            cffT=zeros(2000,2*(siz(2)^2));
        else
            cffT=zeros(2000,(siz(2)^2));
        end
        CxlocT=zeros(500,1);
        CylocT=zeros(500,1);
        
        %go through putative cells, pull pixel region, analyze them
        for j=1:length(maxx);
                yloc=maxy(j);
                xloc=maxx(j);
                if xloc+(siz(1)-1)/2<512+1 & xloc-(siz(1)-1)>0 & yloc+(siz(1)-1)/2<512+1 & yloc-(siz(1)-1)>0 
                   imb=GetBlock(im,im2,siz,xloc,yloc);
                   if std(imb(1:siz(1)^2))/mean(imb(1:siz(1)^2))>0.2
                    [mvt,xy]=lcrosscorr(reshape(imb(1:siz(1)^2),siz(1),siz(1)),circ);
                    yloc=yloc+(mod(xy-1,siz(1))+1)-(siz(1)-1)/2;
                    xloc=xloc+floor(xy/siz(1))+1-(siz(1)-1)/2;
                        if xloc+(siz(1)-1)/2<512+1 & xloc-(siz(1)-1)>0 & yloc+(siz(1)-1)/2<512+1 & yloc-(siz(1)-1)>0 
                            qr=qr+1;
                            mv(qr)=mvt;
                            cffT(qr,:)=GetBlock(im,im2,siz(2),xloc,yloc);
                            CxlocT(qr)=xloc;
                            CylocT(qr)=yloc;
                        end
                   end
                end
        end   
        
        cffT=cffT(1:qr,:);
        CxlocT=CxlocT(1:qr);
        CylocT=CylocT(1:qr);
        
        
        %store the current data, go on to the next loop.
        if storim
            cff(kr+1:kr+qr,:)=cffT;
        else
            for j=1:qr;
                [nf,nmi]=NEnrich(reshape(cffT(j,1:(siz(2)^2)),siz(2),siz(2)),circ);
                cff(kr+j,1:2)=[nf,nmi];
                if yinc; [nf,nmi]=NEnrich(reshape(cffT(j,(siz(2)^2)+1:2*(siz(2)^2)),siz(2),siz(2)),circ); cff(kr+j,3:4)=[nf,nmi]; end
            end    
        end
        Cxloc(kr+1:kr+qr)=CxlocT;
        Cyloc(kr+1:kr+qr)=CylocT;
        tim(kr+1:kr+qr)=i;
        kr=kr+qr;        
    end
    
    %remove matrix entry which were not filled
    cff=cff(1:kr,:);

end

function [nf,nmi]=NEnrich(im,ci)
%compute nuclear enrichment
    siz=size(im,1);
    c2=ci;
    a=conv2(im,c2,'same');
    
    [in,xy]=max(a(:));
    yloc=(mod(xy-1,siz(1))+1);
    xloc=floor(xy/siz(1))+1;
    
    for i=1:siz;
        matx(i,:)=1:siz;
        maty(:,i)=1:siz;
    end
    
    dis=sqrt((xloc-matx(:)).^2+(yloc-maty(:)).^2);
    ne=sort(im(dis<6),'descend');
    nf=mean(ne(1:5))./median(ne);
    nmi=median(ne);
end
    
function [bl]=GetBlock(im,im2,siz,locX,locY)
    %load a portion of an image
    hd=(siz-1)/2;
    tim=im(locX-hd:locX+hd,locY-hd:locY+hd);
    
    if ~isempty(im2)
        rim=im2(locX-hd:locX+hd,locY-hd:locY+hd);
        bl=[tim(:)',rim(:)'];
    else
        bl=tim(:)';
    end

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


function [mv,xy]=lcrosscorr(im,ci)
    siz=size(im,1);


    mv=1;
    %[in,xy]=max(a(:));
    
    imc=conv2(im,ones(5),'same')';
        
    [a,xy]=max(imc(:));

    
end
