function [dat,acqdat]=blockscorrstack(di,p,nframe,siz,circ)

if nargin<5;
    S=load('circ.mat')
    circ=S.circ;
end
%read the acqdat data and strains.txt txt file
[acqdat,strain,numw,numc]=readAcqDat([di,'\acqdat.txt'],[di,'\strains.txt']);


files = dir(strcat(di,'*.tiff'));
asiz=siz^2;
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
    vin=in(find(np(in)==p & (ch(in)=='Y')));
    in=in(find(np(in)==p & (ch(in)=='R')));
    hsiz=floor(siz/2);

    m=0;
    
    imbg=double(imread('bg.png'));
    
    for i=1:nframe

        im=double(imread([di,'\',char(files(in(i)).name)]));
        im2=double(imread([di,'\',char(files(vin(i)).name)]));
        
        %im=im./imbg;
        %im2=im2./imbg;
        
        [cf,x,y]=ublocks(im,im2,siz);
        
%         if i>1
%             for k=1:min(50,length(x))
%                 j=randsample(length(x),1);
%             
%                 a=xcorr2(im(x(j)-hsiz+3:x(j)+hsiz+3,y(j)-hsiz:y(j)+hsiz),imp(x(j)-hsiz:x(j)+hsiz,y(j)-hsiz:y(j)+hsiz)); 
%                 [mv,ml(k)]=max(a(:));
%             
%             end
%             xsh=mean(mod(ml,size(a,1)));
%             ysh=mean(floor(ml./size(a,1))+1);
%             sh(i)=xsh;
%         end
            
        
        %[IDX, K] =
        %kmeans(cf,floor((size(cf,2)/nframes)^(1/1.5)),'distance','correlation');
        
        mv=zeros(1,size(x,2));
        mv2=zeros(1,size(x,2));
        
        for j=1:size(x,2);
            [mv(j),mvin(j)]=lcrosscorr(reshape(cf(j,1:siz^2),[siz,siz]),circ);
            [mv2(j),mvin2(j)]=lcrosscorr(reshape(cf(j,siz^2+1:2*siz^2),[siz,siz]),circ);
        end
        %mvd=sqrt(mod(mvin-1,floor(siz/2)).^2 +  floor(mvin/(floor(siz/2))).^2);
        scor=mv+mv2;
        [k,n]=sort(scor,'descend');
        nmax=floor(length(n)/8);

        x=x(n(1:nmax));
        y=y(n(1:nmax));
        rat=zeros(1,nmax);
        rat2=zeros(1,nmax);
        I=zeros(1,nmax);
        I2=zeros(1,nmax);
                
        for j=1:nmax
                z=sort(cf(n(j),1:siz^2),'descend');
                I(j)=z(floor((siz^2)/4));
                rat(j)=mean(z(1:10))/I(j); 
                z=sort(cf(n(j),siz^2+1:2*siz^2),'descend');
                I2(j)=z(floor((siz^2)/4));
                rat2(j)=mean(z(1:10))/I2(j);     
                %[xc yc sigma] = radialcenter(reshape(cf(n(j),1:asiz),[siz,siz]));
                %x(j)=x(j);
                %y(j)=y(j);
        end
        
        
        dat(1,m+1:m+length(rat))=i;
        dat(2,m+1:m+length(rat))=x;
        dat(3,m+1:m+length(rat))=y;
        dat(4,m+1:m+length(rat))=I;
        dat(5,m+1:m+length(rat))=rat;
        dat(6,m+1:m+length(rat))=I2;
        dat(7,m+1:m+length(rat))=rat2;
        
        m=m+length(rat);
        
        xp=x; yp=y;
        imp=im; imp2=im2;
        
    end
end


function [cim,x,y]=ublocks(im,im2,siz)
    dat=sort(im(:),'descend');
    cutoff=dat(floor(length(im(:))*0.05));

    nbx=floor(size(im,1)/siz);
    nby=floor(size(im,2)/siz);
    cim=zeros(nbx*nby*4,2*siz^2);
    asiz=siz^2;

    j=0;
    for i=0:nbx*nby-1
        tim=im(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        rim=im2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff
            k=k+1;
            if k==length(z); break; end
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz; y(j+1)=floor(i/nby)*siz;
        j=j+1;
    end

    imc=im(floor(siz*0.5):size(im,1),1:size(im,2));
    imc2=im2(floor(siz*0.5):size(im,1),1:size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    for i=0:nbx*nby-1
        tim=imc(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        rim=imc2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff
            k=k+1;
            if k==length(z); break; end
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz+ceil(siz/2); y(j+1)=floor(i/nbx)*siz;
        j=j+1;
    end

    imc=im(1:size(im,1),floor(siz*0.5):size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    imc2=im2(1:size(im,1),floor(siz*0.5):size(im,2));
    for i=0:nbx*nby-1
        tim=imc(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        rim=imc2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff
            k=k+1;
            if k==length(z); break; end
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz; y(j+1)=floor(i/nbx)*siz+ceil(siz/2);
        j=j+1;
    end

    imc=im(floor(siz*0.5):size(im,1),floor(siz*0.5):size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    imc2=im2(floor(siz*0.5):size(im,1),floor(siz*0.5):size(im,2));
    for i=0:nbx*nby-1
        tim=imc(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        rim=imc2(((mod(i,nbx))*siz+1):((mod(i,nbx)+1)*siz),((floor(i/nbx))*siz+1):(floor(i/nbx)+1)*siz);
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff
            k=k+1;
            if k==length(z); break; end
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz+ceil(siz/2); y(j+1)=floor(i/nbx)*siz+ceil(siz/2);
        j=j+1;
    end

    x=x+1+floor(siz/2);
    y=y+1+floor(siz/2);
    np=np./asiz;
    a=find(np>0.2);

    cim=cim(a,:);
    x=x(a);
    y=y(a);
    np=np(a);


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

function [mv,in]=lcrosscorr(im,ci)
    siz=size(im,1);
    %im=conv2(im,ones(3,3),'same');    
    c2=ci(1:13,2:14);
    a=xcorr2(im,c2);
    a2=xcorr2(c2);
    a1=xcorr2(im);
    
    in=max(a(:))./sqrt(max(a1(:))*max(a2(:)));
    mv=in;
%     
%     for i=0:2:floor(siz/2)
%         for j=0:2:floor(siz/2)
%             ai=ci(1:siz-i,1:siz-j);
%             ac=im(i+1:siz,j+1:siz);
%             v(i+1,j+1)=corr(ai(:),ac(:));
%             %v(i+1,j+1)=mean(((ai(:)-mean(ai(:))).*(ac(:)-mean(ac(:)))))/((std(ai(:))*std(ac(:))));
%         end
%     end
%     [mv]=max(v(:));
    
    
end
