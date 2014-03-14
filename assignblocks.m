function [rat,rat2,I,I2]=assignblocks(im,im2, siz,ideal,K)

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
% 
[IDX, C] = kmeans(cim,12,'distance','correlation');
figure
f=find(IDX==1);
length(f)

 for i=1:length(f)
 k=zeros(siz,siz);
 k=reshape(cim(f(i),1:siz^2),siz,siz);
 subplot(ceil(sqrt(length(f))),ceil(sqrt(length(f))),i)
 imagesc(k)
 end


sc=1;
for i=1:length(a)
    sc(i)=corr(cim(i,(1:siz^2))',ideal(1:siz^2)','type','spearman');
end
[z,in]=sort(sc,'descend');

for i=1:50;
    z=sort(cim(in(i),1:siz^2),'descend');
    I(i)=z(floor((siz^2)/2));
    rat(i)=mean(z(1:10))/z(floor((siz^2)/2));
    z=sort(cim(in(i),siz^2+1:2*siz^2),'descend');
    rat2(i)=mean(z(1:10))/z(floor((siz^2)/2));
    I2(i)=z(floor((siz^2)/2));
end


% 
% [z,in]=sort(sc,'descend');
% figure
% for i=1:81
% subplot(9,9,i)
% imagesc(reshape(cim(in(i),1:siz^2),siz,siz))
% end


end
    