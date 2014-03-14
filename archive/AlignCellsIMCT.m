function [cellV,cellx,celly,cfo]=AlignCellsIMCT(cff,Cxloc,Cyloc,tim,ci,nuc)
%function takes output from CovMapCells and aligns cells into traces
%cff is an array of cell images (or NE values)
%Cxloc,Cyloc are arrays of doubles giving X/Y locations of each cell
%tim is an array of doubles giving time stamps to each cell
%ci is a 2s double array giving an ideal cell image
%nuc is a FLAG describing wether to return nuclear enrichment or raw
%intensity of each cell

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

%biuld arrays to store data
if yinc
    cellV=nan(500,2,max(tim));
    if cfi; cfo=nan(5000,2*siz.^2); end
else
    cellV=nan(500,1,max(tim));
    if cfi; cfo=nan(5000,siz.^2); end
end    
    

af=find(tim==min(tim)); mnc=length(af);
[t,tin]=sort(tim);
xs=Cxloc(tin); ys=Cyloc(tin); i=0; na=length(tim);
nsv=max(tim)-min(tim);

%loop while there are cells left to assign
while na>nsv
    i=i+1;
    %pick the first cell in the array as the current cell of interest
    x=Cxloc(1); y=Cyloc(1);
    md=nan(max(tim)-min(tim),1);
    ts=(min(tim):max(tim));
    k=0;
    for j=min(tim):max(tim)
        %step through time and see if you can find cells that are
        %consistent with its localization
        k=k+1;
        a=find(tim==j); mnc=length(a);
        dx=x-Cxloc(a); dy=y-Cyloc(a);
        dt=sqrt(dx.^2+dy.^2);
        mi=find(dt<8); %if at least one cell is within 8 pixels of its previous location
            if length(mi)==1
                mif=mi;
            elseif length(mi)>1
                [mv,mii]=min(dt);
                mif=mii(1);
            else
                mif=[];
            end
            
            if ~isempty(mif)
            %compute nuclear enrichment (or mean intensity), assign to cell
            %data
                [nf,nmi]=NEnrich(reshape(cff(a(mif),1:siz^2),siz,siz),ci);
                if yinc; [nfy,nmiy]=NEnrich(reshape(cff(a(mif),1+siz^2:2*(siz^2)),siz,siz),ci); end
                cellx(i,j)=Cxloc(a(mif)); x=(x+Cxloc(a(mif)))/2;
                celly(i,j)=Cyloc(a(mif)); y=(y+Cyloc(a(mif)))/2;

                %store NE/intensity
                if nuc
                    cellV(i,1,j)=nf; 
                    if yinc; cellV(i,2,j)=nfy; end
                else
                    cellV(i,1,j)=nmi;
                    if yinc; cellV(i,2,j)=nmiy; end
                end
                md(k)=a(mif);
            else
                cellx(i,j)=nan;
                celly(i,j)=nan;
            end
    end
    
    %remove this cell from the dataset to prevent double counting
    md=md(find(md>0));
    a=setxor(1:length(tim),md');
    if ~isempty(a) & min(a)>0
        tim=tim(a);
        Cyloc=Cyloc(a);
        Cxloc=Cxloc(a);
        cff=cff(a,:);
        na=length(tim);
    else
        na=0;
    end

end


cellx=cellx(1:i,:);
celly=celly(1:i,:);
cellV=cellV(1:i,:,:);
mnc=i;

%look through traces, find cells that are capable of being tracked for
%substantial periods of time
for i=1:mnc
    l(i)=sum((cellV(i,1,:))>0);
    maxrl(i)=runlength(cellV(i,1,:));
    maxgl(i)=maxgapsize(cellV(i,1,:));
end
a=find(l>(max(tim)-min(tim))/2 & maxrl>(max(tim)-min(tim))/4);

%return cell data
cellx=cellx(a,:);
celly=celly(a,:);
cellV=cellV(a,:,:);


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
