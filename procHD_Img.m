function [IO,v]=procHD_Img(im,tims)
%im should be a numel(tims),512,512 matrix

for i=1:numel(tims)
    a=im(i,:,:);
    v(i)=quantile(a(:),0.99);
end

%f=fit(1/tims,v,'a*x+b');
scalerat=min(tims)./(tims);

IO=zeros(size(im,2),size(im,3));
j=0;
for i=1:numel(tims)
    a=im(i,:,:);
    if quantile(a(:),0.999)<2^13.8
        IO(:)=IO(:)+a(:)*scalerat(i);
        j=j+1;
    end
end

%IO=IO/j;

end