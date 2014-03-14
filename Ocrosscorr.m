function [mv,xy]=Ocrosscorr(im,ci)
    siz=size(im,1);
    %im=conv2(im,ones(3,3),'same');    
    %c2=ci(1:13,2:14);
    c2=ci;
    a=xcorr2(im,c2);
    a2=xcorr2(c2);
    a1=xcorr2(im);
    
    in=max(a(:))./sqrt(max(a1(:))*max(a2(:)));
    mv=in;
    [in,xy]=max(a(:));
    
    imc=conv2(im,ones(5),'same')';
        
    [a,xy]=max(imc(:));
    
end