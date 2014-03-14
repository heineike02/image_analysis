function [rat,rat2,I,I2]=assayblocks(di,p,nframe,siz,K,poi)

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


        im=imread([di,'\',char(files(in(nframe)).name)]);
        im2=imread([di,'\',char(files(vin(nframe)).name)]);

        cim=ublocks(im,im2,siz);

    try
    [IDX, K] = kmeans(cim,[],'distance','correlation','start',K);

    rat=0; rat2=0; I=0; I2=0;
    k=0;
    for i=1:length(poi)
        a=find(IDX==poi(i));
        if ~isempty(a)
        for j=1:length(a)
            k=k+1;
            z=sort(cim(a(j),1:siz^2),'descend');
            I(k)=z(floor((siz^2)/2));
            rat(k)=mean(z(1:10))/z(floor((siz^2)/2));
            z=sort(cim(a(j),siz^2+1:2*siz^2),'descend');
            rat2(k)=mean(z(1:10))/z(floor((siz^2)/2));
            I2(k)=z(floor((siz^2)/2));
        end
        end
    end
    
    I=mean(I); I2=mean(I2);
    rat=mean(rat); rat2=mean(rat2);
    catch
       I=nan; I2=nan;
       rat=nan; rat2=nan;
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