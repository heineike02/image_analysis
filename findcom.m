function [x,y,s]=findcom(im)

for i=1:size(im,1)
    for j=1:size(im,2)
        bm(1,(i-1)*size(im,1)+j)=i;
        bm(2,(i-1)*size(im,1)+j)=j;
    end
end

rb=[(bm)',(bm.^2)',(bm.^3)',ones(numel(im),1)];

[B,bi,r,rint,stats]=regress(im(:),rb);
s=stats(1);
[a,m]=max(rb*B);

y=mod(m-1,size(im,1))+1;
x=floor(m/size(im,1))+1;