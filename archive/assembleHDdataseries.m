function vs=assembleHDdataseries(di,p,frameHD,tims,siz,nframe,outdi,mchan)

if nargin<8
    mchan='R'; 
end

mkdir(outdi)

%read the acqdat data and strains.txt txt file
[acqdat,strain,numw,numc]=readAcqDat([di,'\acqdat.txt'],[di,'\strains.txt']);

files = dir(strcat(di,'*.tiff')); if isempty(files); files = dir(strcat(di,'\*.tiff')); end
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
    in=in(find(np(in)==p & (ch(in)==mchan(1))));
    
    imbg=double(imread('C:\JSO\MicroscopyData\Background\RFP_BG.tif'));
    
    
    if length(nframe)==1
        framestoload=1:nframe;
    else
        framestoload=nframe;
    end
    
    im=zeros(frameHD,siz,siz);
    
    for i=framestoload
        cf=(i-1)*frameHD+1;
        for j=1:frameHD
            im(j,:,:)=double(imread([di,'\',char(files(in(cf+j-1)).name)]));
        end
        [IO,v]=(procHD_Img(im,tims));
        IO=uint16(IO);
        a=uint16(zeros(512,512));
        a(:,:)=im(1,:,:);
        vs(i,:)=v;
        imwrite(IO,[outdi,'\',mchan(1),'FP_p',num2str(p),'_t',num2str(i),'.tif'],'tiff')    
        imwrite(a,[outdi,'\','QFP_p',num2str(p),'_t',num2str(i),'.tif'],'tiff')    
    end
    
    
        
        