function [cellV,cellx,celly,cfo]=AlignCellsIM(cff,Cxloc,Cyloc,tim,ci,nuc)

cellx=nan(max(tim),5000);
celly=nan(max(tim),5000);
cfi=1;
mitim=min(tim);
matim=max(tim);

if nargin<7; nuc=1; end
if nargin<6; cfi=0; end
if nargin<5 | isempty(ci);
    S=load('circ.mat')
    circ=S.circ;
    ci=circ(1:13,2:14);
end

%test for squareness of input
if mod(sqrt(size(cff,2)),1)~=0
    yinc=1;
    siz=(sqrt(size(cff,2)/2));
else
    yinc=0;
    siz=(sqrt(size(cff,2)));
end

if size(cff,2)==2 | size(cff,2)==4
   storim=0; 
   if size(cff,2)==2
       yinc=0;
   else
       yinc=1;
   end
else
   storim=1;
end


if yinc
    cellV=nan(max(tim),2,500);
    if cfi; cfo=nan(5000,2*siz.^2); end
else
    cellV=nan(max(tim),1,500);
    if cfi; cfo=nan(5000,siz.^2); end
end    
    
    
a=find(tim==min(tim)); mnc=length(a);
for j=1:length(a);
    if storim
        [nf,nmi]=NEnrich(reshape(cff(a(j),1:siz^2),siz,siz),ci);
        if yinc; [nfy,nmiy]=NEnrich(reshape(cff(a(j),1+siz^2:2*(siz^2)),siz,siz),ci); end
    else
        nf=cff(a(j),1);
        nmi=cff(a(j),2);
        if yinc; 
            nfy=cff(a(j),3);
            nmiy=cff(a(j),4); 
        end        
    end
    cellx(min(tim),j)=Cxloc(a(j));
    celly(min(tim),j)=Cyloc(a(j));
    if nuc
        cellV(min(tim),1,j)=nf; 
        if yinc; cellV(min(tim),2,j)=nfy; end
    else
        cellV(min(tim),1,j)=nmi;
        if yinc; cellV(min(tim),2,j)=nmiy; end
    end
    
    if cfi; cfo(j,:)=(cff(a(j),:)); end
end

    for j=1:mnc
        cavX(j)=nanmean(cellx(:,j));
        cavY(j)=nanmean(celly(:,j));
    end
tn=0; tng=10; tne=0;
for i=min(tim)+1:max(tim)
    a=find(tim==i);
    for j=1:length(a);
        if storim
            [nf,nmi]=NEnrich(reshape(cff(a(j),1:siz^2),siz,siz),ci);
            if yinc; [nfy,nmiy]=NEnrich(reshape(cff(a(j),1+siz^2:2*(siz^2)),siz,siz),ci); end
        else
            nf=cff(a(j),1);
            nmi=cff(a(j),2);
            if yinc; 
                nfy=cff(a(j),3);
                nmiy=cff(a(j),4); 
            end  
        end
            dx=cellx(i-1,1:mnc)-Cxloc(a(j));
            dy=celly(i-1,1:mnc)-Cyloc(a(j));

            dt=sqrt(dx.^2+dy.^2);
            mi=find(dt<9);
            if length(mi)==1
                mif=mi;
                tn=tn+1;
            elseif length(mi)>1
                dz=0;
                for k=1:length(mi); dz(k)=max(max(conv2(reshape(cfo(mi(k),:),17,17*2),reshape(cff(a(j),:),17,17*2)))); end
                dz=max(dz)-dz;
                dt=sqrt(dx(mi).^2+dy(mi).^2+dz.^2);
                [mv,mii]=min(dt);
                mif=mi(mii);
            else
                dx=cavX-Cxloc(a(j));
                dy=cavY-Cyloc(a(j));
                dt=sqrt(dx.^2+dy.^2);
                [mv,mi]=min(dt);
                if mv<12
                    mi=floor(mi./(i-1))+1;
                    mif=mi;
                    tng=tng+1;
                    mif=[];
                else
                    mnc=mnc+1;
                    mif=mnc;
                    tne=tne+1;
                end
            end
            if ~isempty(mif)
                cellx(i,mif)=Cxloc(a(j));
                celly(i,mif)=Cyloc(a(j));

                if nuc
                    cellV(i,1,mif)=nf; 
                    if yinc; cellV(i,2,mif)=nfy; end
                else
                    cellV(i,1,mif)=nmi;
                    if yinc; cellV(i,2,mif)=nmiy; end
                end

                if cfi; cfo(mif,:)=(cff(a(j),:)); end
            end
        
    end
    %nanmean(cellx(i,:)-cellx(i-1,:))
    %nanmean(celly(i,:)-celly(i-1,:))
    for j=1:mnc
        cavX(j)=nanmean(cellx(mitim:i,j));
        cavY(j)=nanmean(celly(mitim:i,j));
        if isempty(cellx(i,j))
            cellx(i,j)=cellx(i-1,j);
            celly(i,j)=celly(i-1,j);
        end
    end
end

cellx=cellx(:,1:mnc);
celly=celly(:,1:mnc);
cellV=cellV(:,:,1:mnc);
if cfi; cfo=cfo(1:mnc,:); end

for i=1:mnc
    l(i)=sum(~isnan(cellV(:,1,i)));
    maxrl(i)=runlength(cellV(:,1,i));
    maxgl(i)=maxgapsize(cellV(:,1,i));
end

a=find(l>(max(tim)-min(tim))/2 & maxrl>(max(tim)-min(tim))/4);


cellx=cellx(:,a);
celly=celly(:,a);
cellV=cellV(:,:,a);
if cfi;  cfo=cfo(:,a,:); else; cfo=[]; end


end


function maxrl=runlength(dat)
maxrl=0;    
k=0;
smrl=0;
for i=1:length(dat)
        if ~isnan(dat(i))
            maxrl=maxrl+1;
        else
            if maxrl>0
                k=k+1;
                smrl(k)=maxrl;
            end
            maxrl=0;
        end
end
k=k+1;
smrl(k)=maxrl;

maxrl=max(smrl);
end

function maxrl=maxgapsize(dat)
maxrl=0;    
k=0;
smrl=0;
for i=1:length(dat)
        if isnan(dat(i))
            maxrl=maxrl+1;
        else
            if maxrl>0
                k=k+1;
                smrl(k)=maxrl;
            end
            maxrl=0;
        end
end
k=k+1;
smrl(k)=maxrl;

maxrl=max(smrl);
end


function [nf,nmi]=NEnrich(im,ci)
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
    nmi=median(ne);
end

function m=fastnanmean(x)
    m=mean(x(~isnan(x)));
end
