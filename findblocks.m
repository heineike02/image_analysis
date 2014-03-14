function [K,poi]=findblocks(di,p,nframes,siz,circ)

addpath('C:\program files\MATLAB\r2009a\toolbox\stats')

files = dir(strcat(di,'*.tiff'));

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



    for i=1:nframes
        im=imread([di,'\',char(files(in(i)).name)]);
        im2=imread([di,'\',char(files(vin(i)).name)]);

        cim=ublocks(im,im2,siz);
        
        if i==1
            cf=cim';
        else
            cf=[cf,cim'];
        end
    end

    cf=cf';
    
    [IDX, K] = kmeans(cf,floor((size(cf,2)/nframes)^(1/1.5)),'distance','correlation');
    
    for i=1:size(cf,1);
        a=(xcorr2(circ,reshape(cf(i,1:17^2),[17,17]))); [sm(i),l(i)]=max(a(:)); l(i)=sqrt((16-mod(l(i),33))^2+(16-l(i)/33)^2); s(i)=sum(a(:));
    end
    
    
    if nargin==5
        for i=1:size(K,1);
            a=(xcorr2(circ,reshape(K(i,1:17^2),[17,17]))); [sm(i),l(i)]=max(a(:)); l(i)=sqrt((16-mod(l(i),33))^2+(16-l(i)/33)^2); s(i)=sum(a(:));
        end

        poi=find(l<3 & sm>.3);
    else
        poi=[];
    end
    
    
end



function cim=ublocks(im,im2,siz)
    dat=sort(im(:),'descend');
    cutoff=dat(floor(length(im(:))*0.05));

    nbx=floor(size(im,1)/siz);
    nby=floor(size(im,2)/siz);
    cim=zeros(nbx*nby*4,2*siz^2);

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
        j=j+1;
    end

    a=find(np>0.2*(siz^2));

    cim=cim(a,:);


end