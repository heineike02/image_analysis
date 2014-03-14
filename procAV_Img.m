function [IO]=procAV_Img(im,tims)
%im should be a numel(tims),512,512 matrix

IO=zeros(size(im,2),size(im,3));
j=0;
for i=1:numel(tims)
    a=im(i,:,:);
    IO(:)=IO(:)+a(:)*tims(i);
end

%IO=IO/j;

end