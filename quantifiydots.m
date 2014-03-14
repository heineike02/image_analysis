function imax=quantifiydots(di)

files = dir(strcat(di,'*.tif'));

im=zeros(length(files),512,512);


for i=1:length(files);
   im(i,:,:)=(imread([di,'\',char(files(i).name)]));
end

imm=max(im);

[v,li]=sort(imm(:),'descend');



for i=1:100;
    lx=mod((li(i)-1),512);
    ly=(floor((li(i)-1)/512))+1;
    rx=min(512,lx+4); lx=max(1,lx-4); 
    ry=min(512,ly+4); ly=max(1,ly-4); 

    
    la=zeros(rx-lx+1,ry-ly+1);
    la(:,:)=imm(1,lx:rx,ly:ry);
    
    imax(i)=max(la(:));
    
    imm(1,lx:rx,ly:ry)=0;
    [v,li]=sort(imm(:),'descend');
end


end