function NE=ProcessCells(cff,ci)


%test for squareness of input
if mod(sqrt(size(cff,2)),1)~=0
    yinc=1;
    siz=(sqrt(size(cff,2)/2));
else
    yinc=0;
    siz=(sqrt(size(cff,2)));
end


for i=1:size(cff,1);
    NE(i,1)=NEnrich(reshape(cff(i,1:siz^2),siz,siz),ci);
    if yinc; NE(i,2)=NEnrich(reshape(cff(i,1+siz^2:2*(17^2)),siz,siz),ci); end
end

end


function nf=NEnrich(im,ci)
    siz=size(im,1);
    %im=conv2(im,ones(3,3),'same');    
    %c2=ci(1:13,2:14);
    c2=ci;
    a=conv2(im,c2,'same');
    
    [in,xy]=max(a(:));
    yloc=(mod(xy-1,siz(1))+1);
    xloc=floor(xy/siz(1))+1;
    
    for i=1:siz;
        matx(i,:)=1:siz;
        maty(:,i)=1:siz;
    end
    
    dis=sqrt((xloc-matx(:)).^2+(yloc-maty(:)).^2);
    ne=sort(im(dis<6),'descend');
    nf=mean(ne(1:5))./median(ne);
    
end