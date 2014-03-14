function [cff,acqdat,Cxloc,Cyloc,tim]=FindUBlocks(di,p,nframe,siz,circ)

if nargin<5;
    S=load('circ.mat')
    circ=S.circ;
    circ=circ(1:13,2:14);
end
%read the acqdat data and strains.txt txt file
[acqdat,strain,numw,numc]=readAcqDat([di,'\acqdat.txt'],[di,'\strains.txt']);


files = dir(strcat(di,'*.tiff'));
asiz=siz(1)^2;
k=0; kr=0;

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
    vin=in(find(np(in)==p & (ch(in)=='Y')));
    in=in(find(np(in)==p & (ch(in)=='R')));
    hsiz=floor(siz(1)/2);

    m=0;
    
    if isempty(vin); yinc=0;else yinc=1; end

    m=0;
    
    imbg=double(imread('C:\JSO\MicroscopyData\Background\RFP_BG.tif'));
    kr=0; qr=0;
    
    if yinc
        cff=zeros(10000,2*(siz(2)^2));
    else
        cff=zeros(10000,(siz(2)^2));
    end
    
        
    if length(nframe)==1
        framestoload=1:nframe;
    else
        framestoload=nframe;
    end
    
    
    for i=framestoload

        im=double(imread([di,'\',char(files(in(i)).name)]));
        
        
        if yinc
            im2=double(imread([di,'\',char(files(vin(i)).name)]));
            im2=im2./imbg;
        else
            im2=[];
        end
        
        im=im./imbg;
        
        [cf,x,y]=ublocks(im,im2,siz(1));
        

        %[IDX, K] =
        %kmeans(cf,floor((size(cf,2)/nframes)^(1/1.5)),'distance','correlation');
        
        mv=zeros(1,size(x,2));
        mv2=zeros(1,size(x,2));
        
        for j=1:size(x,2);
            [mv(j),mvin(j)]=lcrosscorr(reshape(cf(j,1:siz(1)^2),[siz(1),siz(1)]),circ);
            if yinc; [mv2(j),mvin2(j)]=lcrosscorr(reshape(cf(j,siz(1)^2+1:2*siz(1)^2),[siz(1),siz(1)]),circ); end
        end
        %mvd=sqrt(mod(mvin-1,floor(siz/2)).^2 +  floor(mvin/(floor(siz/2))).^2);
        if yinc;
            scor=mv+mv2;
        else
            scor=mv;
        end
        [k,n]=sort(scor,'descend');
        nmax=floor(length(n));
                
       
         qr=0;
        if yinc
            cffT=zeros(200,2*(siz(2)^2));
        else
            cffT=zeros(200,(siz(2)^2));
        end
        CxlocT=zeros(500,1);
        CylocT=zeros(500,1);
        mv=zeros(500,1);
        for j=1:nmax;
                yloc=y(n(j))+(mod(mvin(n(j))-1,siz(1))+1)-(siz(1)-1)/2;
                xloc=x(n(j))+floor(mvin(n(j))/siz(1))+1-(siz(1)-1)/2;
                if xloc+(siz(1)-1)/2<512+1 & xloc-(siz(1)-1)>0 & yloc+(siz(1)-1)/2<512+1 & yloc-(siz(1)-1)>0 
                   imb=GetBlock(im,im2,siz(1),xloc,yloc);
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
        mv=mv(1:qr);
        
        dis=squareform(pdist([CxlocT,CylocT]));
        dis(eye(length(dis))>0)=max(dis(:));
        
        k=0;
        while min(dis(:))<5
            dnu=find(dis(k+1,:)>5);
            CxlocT=CxlocT(dnu);
            CylocT=CylocT(dnu);           
            cffT=cffT(dnu,:);
            mv=mv(dnu);                  
            k=k+1;
            qr=length(dnu);
            dis=squareform(pdist([CxlocT,CylocT]));
            dis(eye(length(dis))>0)=max(dis(:));
        end
            
        
        cff(kr+1:kr+qr,:)=cffT;
        Cxloc(kr+1:kr+qr)=CxlocT;
        Cyloc(kr+1:kr+qr)=CylocT;
        tim(kr+1:kr+qr)=i;
        kr=kr+qr;
        
        
    end
    
    cff=cff(1:kr,:);
    
end

function [bl]=GetBlock(im,im2,siz,locX,locY)

    hd=(siz-1)/2;
    tim=im(locX-hd:locX+hd,locY-hd:locY+hd);
    
    if ~isempty(im2)
        rim=im2(locX-hd:locX+hd,locY-hd:locY+hd);
        bl=[tim(:)',rim(:)'];
    else
        bl=tim(:)';
    end

end


function [cim,x,y]=ublocks(im,im2,siz)
    dat=sort(im(:),'descend');
    cutoff=dat(floor(length(im(:))*0.2));

    nbx=floor(size(im,1)/siz);
    nby=floor(size(im,2)/siz);
    if ~isempty(im2)
        cim=zeros(nbx*nby*4,2*siz^2);
    else
        cim=zeros(nbx*nby*4,siz^2);
    end
        asiz=siz^2;

    j=0;
    for i=0:nbx*nby-1
        tim=im(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        if ~isempty(im2)
            rim=im2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
            cim(j+1,:)=[tim(:)',rim(:)'];
        else
            cim(j+1,:)=[tim(:)'];
        end
        z=sort(tim(:),'descend');
        k=0;
        nps(j+1)=std(tim(:))./mean(tim(:));
        x(j+1)=mod(i,nbx)*siz+1; y(j+1)=floor(i/nby)*siz+1;
        j=j+1;
    end

    imc=im(floor(siz*0.5):size(im,1),1:size(im,2));
    if ~isempty(im2); imc2=im2(floor(siz*0.5):size(im,1),1:size(im,2)); end
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    for i=0:nbx*nby-1
        tim=imc(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        if ~isempty(im2)
        rim=imc2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
            cim(j+1,:)=[tim(:)',rim(:)'];
        else
            cim(j+1,:)=[tim(:)'];
        end
        z=sort(tim(:),'descend');
        k=0;
        nps(j+1)=std(tim(:))./mean(tim(:));
        x(j+1)=mod(i,nbx)*siz+ceil(siz/2)+1; y(j+1)=floor(i/nbx)*siz+1;
        j=j+1;
    end

    imc=im(1:size(im,1),floor(siz*0.5):size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    if ~isempty(im2); imc2=im2(1:size(im,1),floor(siz*0.5):size(im,2)); end
    for i=0:nbx*nby-1
        tim=imc(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        if ~isempty(im2)
            rim=imc2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
            cim(j+1,:)=[tim(:)',rim(:)'];
        else
            cim(j+1,:)=[tim(:)'];
        end
        z=sort(tim(:),'descend');
        k=0;
        nps(j+1)=std(tim(:))./mean(tim(:));
        x(j+1)=mod(i,nbx)*siz+1; y(j+1)=floor(i/nbx)*siz+ceil(siz/2)+1;
        j=j+1;
    end

    imc=im(floor(siz*0.5):size(im,1),floor(siz*0.5):size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    if ~isempty(im2); imc2=im2(floor(siz*0.5):size(im,1),floor(siz*0.5):size(im,2)); end
    for i=0:nbx*nby-1
        tim=imc(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        if ~isempty(im2)
            rim=imc2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
            cim(j+1,:)=[tim(:)',rim(:)'];
        else
            cim(j+1,:)=[tim(:)'];
        end
        z=sort(tim(:),'descend');
        k=0;
        nps(j+1)=std(tim(:))./mean(tim(:));
        x(j+1)=mod(i,nbx)*siz+ceil(siz/2)+1; y(j+1)=floor(i/nbx)*siz+ceil(siz/2)+1;
        j=j+1;
    end

    x=x+1+floor(siz/2);
    y=y+1+floor(siz/2);
    a=find(nps>0.20);

    cim=cim(a,:);
    x=x(a);
    y=y(a);


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


function AlignImages(xp,yp,x,y)

    
    for i=1:7;
        iold=zeros(512,512);
        iold(xp,yp)=1;
        inew=zeros(512,512);
        yc=y-36+i*9; c=find(yc>0);
        xc=x(c); yc=yc(c);
        inew(xc,yc)=1;
        z(i,:)=xcorr(iold(:),inew(:),27);
    end
    
    [a,l]=max(z(:));
    lou.x=mod(l,55);
    lou.y=1+floor(l/55);
    
end

function [mv,xy]=lcrosscorr(im,ci)
    siz=size(im,1);
    %im=conv2(im,ones(3,3),'same');    
    %c2=ci(1:13,2:14);
    c2=ci;
    a=xcorr2(im,c2);
    a2=xcorr2(c2);
    a1=xcorr2(im);
    
    in=max(a(:))./sqrt(max(a1(:))*max(a2(:)));
    mv=in;
    [in,xy]=max(a(:));
    
    imc=conv2(im,ones(5),'same')';
        
    [a,xy]=max(imc(:));

    
end
