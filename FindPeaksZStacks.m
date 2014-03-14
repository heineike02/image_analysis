function [a]=FindPeaksZStacks(di)

for i=1:5; for j=1:5;
        circ(i,j)=sqrt((i-3)^2+(j-3)^2);
end; end
circ=pdf('norm',circ,0,0.5);

circ=double(imread('C:\JSO\matfiles\work\imblock\deconv.tif'));
circ=circ/max(circ(:));

file = dir(strcat(di,'*.tiff')); file2 = dir(strcat(di,'*.tif'));
files=[file',file2'];
    imbg=double(imread('C:\JSO\MicroscopyData\Background\RFP_BG.tif'));
    imbg=imbg./mean(imbg(:));
    
    im=zeros(length(files),512,512);
    
    for i=1:length(files)
        im(i,:,:)=double(imread([di,'\',char(files(i).name)]));
    end
    
    im=max(im,[],1);
    k=zeros(512,512); k(:,:)=im;
    ima=deconvlucy((k),circ,10);
    %imc=ima./conv2(imbg,ones(9),'same');
    imc=conv2(ima,circ,'same');
    thr=median(imc(:))+4*std(imc(:));

    [px,py]=find(imc>thr);
    imc=ima./conv2(imbg./conv2(imbg,ones(9),'same'),ones(9),'same');
    imc=conv2(imc,ones(5),'same');
       
    [maxx,maxy]=FindMaxima2(imc,px,py,1);
    [maxx,maxy]=FindMaxima2(imc,maxx,maxy,6);
    %[maxx,maxy]=FindLocPoint(ima,maxx,maxy,15,2);
    for i=1:length(maxx); 
        a(i)=sum(GetBlock(k,[],5,maxx(i),maxy(i)));
        %a(i)=k(maxx(i),maxy(i)); 
    end;
    a=a-25*median(k(:));
     figure(1); imagesc(k)
     figure(2); imd=zeros(512); for i=1:length(maxx); imd(maxx(i),maxy(i))=1; end; imagesc(conv2(imd,ones(5)))
     figure(3); hist(a)
    
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