function [mx,my]=AlignImg(di,p,nframe,siz)

mchan='RY';  

file = dir(strcat(di,'*.tiff')); file2 = dir(strcat(di,'*.tif'));
files=[file',file2'];

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
    
    if isempty(vin); yinc=0;else yinc=1; end
    
    if length(nframe)==1
        framestoload=2:nframe;
    else
        framestoload=nframe(2:length(nframe));
    end
    
    mx(1)=0;
    my(1)=0;
    
    for i=framestoload

        im=double(imread([di,'\',char(files(in(i-1)).name)]));
        im2=double(imread([di,'\',char(files(in(i)).name)]));
    
        for j=1:200
            xloc=floor(rand*(512-33))+18;
            yloc=floor(rand*(512-33))+18;
            imb=GetBlock(im,[],31,xloc,yloc);
            for k=1:3; 
                for m=1:3;
                imc=GetBlock(im2,[],31,xloc+k-2,yloc+m-2);
                crc(k,m)=corr(imc',imb');
                end            
            end
            [mv(j),v]=max(crc(:));
            mv(j)=mv(j)./mean(crc(:));
            yoff(j)=floor((v-1)/3)-1;
            xoff(j)=mod(v-1,3)-1;
        end
        
        mx(i)=mean(xoff);
        my(i)=mean(yoff);
        
    end

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


function [mv,xy]=lcrosscorr(im,ci)
    siz=size(im,1);
    mv=1;
    imc=conv2(im,ones(5),'same')';
    [a,xy]=max(imc(:));
    a=a/mean(imc(:));
end