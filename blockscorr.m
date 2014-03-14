function [rat,ratA]=blockscorr(im,im2,siz,circ)
%function [rat,rat2]=blockscorr(di,p,nframe,siz,circ)
%Function blockscorr
%   splits an image into sqares with side length given by 'siz'
%   im/im2 are integer or double images of the same size
%   each square is assayed for containing sufficient numbers of non-zero
%   pixels (use a 20% cutoff)
%   squares that pass this cuttoff are cross correlated with a image given
%   by 'circ' -- typically a single cell
%   squares in the top 12.5% of correlation scores are scored for nuclear
%   localization
    asiz=siz^2;
    k=0;

    for i=1:siz
        matX(i,:)=[1:siz];
        matY(:,i)=[1:siz];
    end
    %split image into blocks, return blocks with objects
    [cf,x,y]=ublocks(im,im2,siz);
    
    %run cross corr on both channels with circ image, find max value
    for i=1:size(cf,1);
        a=(xcorr2(circ,reshape(cf(i,1:asiz),[siz,siz]))); [sm(i),l(i)]=max(a(:)); sm(i)=sm(i)./sum(a(:)); l(i)=sqrt((16-mod(l(i),33))^2+(16-l(i)/33)^2);
        a=(xcorr2(circ,reshape(cf(i,asiz+1:2*asiz),[siz,siz]))); [sm2(i),l2(i)]=max(a(:));   sm2(i)=sm2(i)./sum(a(:)); l2(i)=sqrt((16-mod(l2(i),33))^2+(16-l2(i)/33)^2);
            [xc yc sigma] = radialcenter(reshape(cf(i,1:asiz),[siz,siz]));
            dismat=sqrt((matX-xc).^2+(matY-yc).^2);
            
            [nk,nkp]=sort(dismat(:));

            cm=(cumsum(cf(i,nkp)));
            dq(i)=quantile(cm,0.1);
            pq(i)=quantile(cm,0.9);
            df(i)=sum(cm)/(length(nkp)*sum(cf(i,nkp)));
    end
    
    
    
    %sum correlation of both images to reduce noise
    scor=sm+sm2;
    [k,n]=sort(scor,'descend');
    nmax=floor(length(n));
    

    %nuclear loc score of nmax images
    for i=1:nmax
            z=sort(cf(n(i),1:asiz),'descend');
            rat(i)=mean(z(1:10))/z(floor((asiz)/4)); 
            z=sort(cf(n(i),asiz+1:2*asiz),'descend');
            rat2(i)=mean(z(1:10))/z(floor((asiz)/4));
            
            
            
            [xc yc sigma] = radialcenter(reshape(cf(n(i),1:asiz),[siz,siz]));
            dismat=sqrt((matX-xc).^2+(matY-yc).^2);
            
            [nk,nkp]=sort(dismat(:));

            cm=(cumsum(cf(n(i),nkp)));
            dq(i)=quantile(cm,0.1);
            pq(i)=quantile(cm,0.8);
            
            np=find(dismat<=2); cp=find(dismat>2 & dismat<4);
            ratA(i)=mean(cf(n(i),np))./mean(cf(n(i),cp));
            ratB(i)=mean(cf(n(i),siz^2+np))./mean(cf(n(i),siz^2+cp));
            
    end
    
    figure; plot(dq./pq)
    
    rat=mean(rat); rat2=mean(rat2); ratA=mean(ratA);
end


function [cim,x,y]=ublocks(im,im2,siz)

%function ublocks
%   splits an image into squares with side length 'siz' with four different
%   frames

    dat=sort(im(:),'descend');
    im=im-quantile(dat,0.95);
    
    dat=sort(im(:),'descend');
    cutoff=dat(floor(length(im(:))*0.05));

    asiz=siz^2;
    j=0;

    nbx=floor(size(im,1)/siz);
    nby=floor(size(im,2)/siz);
    cim=zeros(nbx*nby*4,2*asiz);
    
    np=zeros(nbx*nby*4,1);
    x=zeros(nbx*nby*4,1);
    y=zeros(nbx*nby*4,1);
    
    ia=0:nbx*nby-1;
    
    xposl=((mod(ia,nbx))*siz+1);  
    xposh=((mod(ia,nbx)+1)*siz);
    yposl=((floor((ia)/nbx))*siz+1);
    yposh=(floor(ia/nbx)+1)*siz;
    
    for i=1:nbx*nby
        tim=im(xposl(i):xposh(i),yposl(i):yposh(i));
        rim=im2(xposl(i):xposh(i),yposl(i):yposh(i));
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff && k<asiz-2
            k=k+1;
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz+ceil(siz/2); 
        y(j+1)=floor(i/nbx)*siz+ceil(siz/2);
        j=j+1;
    end

    imc=im(floor(siz*0.5):size(im,1),1:size(im,2));
    imc2=im2(floor(siz*0.5):size(im,1),1:size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    
    xposl=((mod(ia,nbx))*siz+1);  
    xposh=((mod(ia,nbx)+1)*siz);
    yposl=((floor((ia)/nbx))*siz+1);
    yposh=(floor(ia/nbx)+1)*siz;
    
    for i=1:nbx*nby
        tim=imc(xposl(i):xposh(i),yposl(i):yposh(i));
        rim=imc2(xposl(i):xposh(i),yposl(i):yposh(i));
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff && k<asiz-2
            k=k+1;
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz+2*ceil(siz/2); y(j+1)=floor(i/nbx)*siz+ceil(siz/2);
        j=j+1;
    end

    imc=im(1:size(im,1),floor(siz*0.5):size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    imc2=im2(1:size(im,1),floor(siz*0.5):size(im,2));
    
    xposl=((mod(ia,nbx))*siz+1);  
    xposh=((mod(ia,nbx)+1)*siz);
    yposl=((floor((ia)/nbx))*siz+1);
    yposh=(floor(ia/nbx)+1)*siz;
    
    for i=1:nbx*nby
        tim=imc(xposl(i):xposh(i),yposl(i):yposh(i));
        rim=imc2(xposl(i):xposh(i),yposl(i):yposh(i));
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff && k<asiz-2
            k=k+1;
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz+ceil(siz/2); y(j+1)=floor(i/nbx)*siz+2*ceil(siz/2);
        j=j+1;
    end

    imc=im(floor(siz*0.5):size(im,1),floor(siz*0.5):size(im,2));
    nbx=floor(size(imc,1)/siz);
    nby=floor(size(imc,2)/siz);
    imc2=im2(floor(siz*0.5):size(im,1),floor(siz*0.5):size(im,2));
    
    xposl=((mod(ia,nbx))*siz+1);  
    xposh=((mod(ia,nbx)+1)*siz);
    yposl=((floor((ia)/nbx))*siz+1);
    yposh=(floor(ia/nbx)+1)*siz;
    
    for i=1:nbx*nby
        tim=imc(xposl(i):xposh(i),yposl(i):yposh(i));
        rim=imc2(xposl(i):xposh(i),yposl(i):yposh(i));
        cim(j+1,:)=[tim(:)',rim(:)'];
        z=sort(tim(:),'descend');
        k=0;
        while z(k+1)>cutoff && k<asiz-2
            k=k+1;
        end
        np(j+1)=k;
        x(j+1)=mod(i,nbx)*siz+2*ceil(siz/2); y(j+1)=floor(i/nbx)*siz+2*ceil(siz/2);
        j=j+1;
    end

    a=find(np>0.2*(asiz));

    cim=cim(a,:);
    x=x(a);
    y=y(a);
end
