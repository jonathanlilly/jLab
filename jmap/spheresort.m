function[varargout]=spheresort(varargin)
%SPHERESORT  Sorted great circle distances to nearby points on the earth.
%
%   Computing great circle distances on the earth between two sets of
%   points --- data points and a fixed grid, say --- is computationally 
%   expensize.  SPHERESORT speeds this up substantially. 
%
%   [DS,XS,YS]=SPHERESORT(LAT,LON,LATO,LONO,CUTOFF) returns the great
%   circle distances DS between data points at locations LAT, LON and 
%   nearby grid points located at LATO, LONO, sorting in order of distance. 
%
%   XS and YS are the corresponding positions of the sorted data points 
%   in a local tangent plane about each grid point. 
%
%   SPHERESORT only computes distances for nearby points.  CUTOFF is the
%   maximum distance (in kilometers) for which we wish the great circle
%   distance to be computed.  
%
%   CUTOFF should be less that 1/4 of the circumference of the earth, 
%   RADEARTH * PI/2, so that (LAT,LON) points lie in the same hemisphere
%   as the corresponding grid points.
% 
%   LAT and LON are arrays of the same size into data point locations.
%   LATO and LONO are arrays of length M and N, say, specifying the
%   latitudes and longitudes of an M x N matrix of grid points, i.e.
%
%       LATO =  [LATO_1;    LONO= [LONO_1 LONO_2 ... LONO_N]. 
%                LATO_2; 
%                  ...
%                LATO_M]
%
%   The output arrays are then each M x N x P arrays of column vectors, 
%   where P is the maximum number of points in the CUTOFF neighborhood at
%   any grid point.  Points farther away are filled with NANs.
%
%   DS gives the distances of data points less than CUTOFF kilometers 
%   from the (m,n)th grid point, sorted in order of increasing distance.  
%
%   XS and YS are also M x N x P, and give the coordinates of the 
%   corresponding data points in a Cartesian expansion about the (m,n)th
%   grid point, in kilometers.  See LATLON2XY for details.
%   _________________________________________________________________
% 
%   Additional input parameters
%
%   Let's say some additional variables Z1, Z2,...,ZK are given at
%   the data locations LAT, LON.  Then 
%   
%   [DS,XS,YS,Z1S,Z2S,...,ZKS]=
%
%                SPHERESORT(LAT,LON,Z1,Z2,...,ZK,LATO,LONO,CUTOFF);
%
%   also returns the sorted values of these variables.
%
%   Z1S, Z2S,...,ZKS are the same size as the other output arguments.  
%   Z1S then gives the value of Z1 at data points no more than CUTOFF 
%   kilometers from the (m,n)th grid point, etc.
%
%   To output an index into the data point locations, use
%  
%   [DS,XS,YS,INDEX]=SPHERESORT(LAT,LON,1:LENGTH(LAT(:)),LATO,LONO,CUTOFF).
%   _________________________________________________________________
%  
%   Parellel algorithm
%
%   With Matlab's Parallel Computing Toolbox installed, SPHERESORT can 
%   take advantage of multiple cores to speed up operations.
%
%   SPHERESORT('parallel',NWORKERS) parallelizes the computation with the
%   number of workers set to NWORKERS.  Choose this as the number of 
%   independent processors on your machine. 
%
%   As the additional efficiency is not dramatic, this is typically only an 
%   advantage for very large datasets.  Due to the nature of the
%   calculation we have to use SPMD, which is less efficient than PARFOR.
%
%   As an example, for a dataset with 1 million points, a 12 core Mac Pro 
%   takes about 44 seconds to complete the sort on a 1x1 degree grid, 
%   versus 184 seconds for the standard algorithm, a factor of 4 faster.
%   _________________________________________________________________
%
%   See also TWODSORT, POLYSMOOTH.
%
%   'spheresort --t' runs a test.
%
%   Usage: ds=spheresort(lat,lon,lato,lono,cutoff);
%          ds=spheresort(lat,lon,lato,lono,cutoff,'parallel',Nworkers);
%          [ds,xs,ys]=spheresort(lat,lon,lato,lono,cutoff);
%          [ds,xs,ys,z1s,z2s]=spheresort(lat,lon,z1,z2,lato,lono,cutoff);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details
 

if strcmpi(varargin{1}, '--t')
    spheresort_test;
    %spheresort_parallel_test;
    return
end

lat=varargin{1};
lon=varargin{2};

algstr='current';
str='serial';
Nworkers=[];

for i=1:2
    if ischar(varargin{end-1})
        str=varargin{end-1};
        Nworkers=varargin{end};
        varargin=varargin(1:end-2);
    elseif ischar(varargin{end})
        algstr=varargin{end};
        varargin=varargin(1:end-1);
    end
end

P=inf;
    
lato=varargin{end-2};
lono=varargin{end-1};
cutoff=varargin{end};
args=varargin(3:end-3);

if ~aresame(size(lat),size(lon))
    error('LAT and LON must be the same size.')
end

cutoff=abs(cutoff);
if cutoff>(radearth*pi/2)
    error('Sorry, SPHERESORT requires CUTOFF<RADEARTH*PI/2.')
end

for i=1:nargin
    varargout{i}=[];
end
indexo=find(isfinite(lat)&isfinite(lon));
if ~isempty(indexo)
    vcolon(lat,lon);
    vindex(lat,lon,indexo,1);
else
    disp(['No finite data values.']), return
end

vcolon(lato,lono);
lono=lono';

if strcmpi(str(1:3),'par')
    [d,xd,yd,indexd]=spheresort_spmd(lat,lon,lato,lono,cutoff,Nworkers);
else
    if strcmpi(algstr(1:3),'cur')
        [d_cell,xd_cell,yd_cell,indexd_cell]=spheresort_current(lat,lon,lato,lono,cutoff);
    else
        [d_cell,xd_cell,yd_cell,indexd_cell]=spheresort_former(lat,lon,lato,lono,cutoff);
    end
    
    disp('SPHERESORT reorganizing, this may take a while...')
    
    L=cellength(d_cell);    
    [d,xd,yd,indexd]=vzeros(size(L,1),size(L,2),maxmax(L),nan);
    for i=1:size(L,1)
        for j=1:size(L,2)
            if L(i,j)>=1
                d(i,j,1:L(i,j))=d_cell{i,j};
                xd(i,j,1:L(i,j))=xd_cell{i,j};
                yd(i,j,1:L(i,j))=yd_cell{i,j};
                indexd(i,j,1:L(i,j))=indexd_cell{i,j};
            end
        end
    end
    clear d_cell xd_cell yd_cell indexd_cell
end

varargout{1}=d;
varargout{2}=xd;
varargout{3}=yd;

%Finally, account for extra input arguments
for k=1:length(args)
    temp=nan*d;
    nonnani=~isnan(indexd);
    temp(nonnani)=args{k}(indexo(indexd(nonnani)));
    varargout{k+3}=temp;
end

disp('SPHERESORT finished.')




function [d,xd,yd,indexd]=spheresort_spmd(lat,lon,lato,lono,cutoff,Nworkers)

%Parallelize with single program multiple data
dlat_cutoff=jrad2deg(cutoff./radearth);
b=floor((length(lato)/Nworkers)*[1:Nworkers]);
a=[1 b(1:end-1)+1];
clear struct labindex
index=[1:length(lat)];
spmd
    %Determine the data to send to each worker
    spmdindex=a(labindex):b(labindex);
    bool1=(lat-max(lato(spmdindex)))<=dlat_cutoff;
    bool2=(min(lato(spmdindex))-lat)<=dlat_cutoff;
    bool=bool1&bool2;
    [dp,xdp,ydp,indexdp]=spheresort(lat(bool),lon(bool),index(bool),lato(spmdindex),lono,cutoff);
end

N=zeros(length(dp),1);
for i=1:length(dp)
    N(i)=size(dp{i},3);
end
N=max(N);

[d,xd,yd,indexd]=vzeros(length(lato),length(lono),N,nan);
n=0;
for i=1:length(dp)
    d(n+(1:size(dp{i},1)),:,1:size(dp{i},3))=dp{i};
    xd(n+(1:size(dp{i},1)),:,1:size(dp{i},3))=xdp{i};
    yd(n+(1:size(dp{i},1)),:,1:size(dp{i},3))=ydp{i};
    indexd(n+(1:size(dp{i},1)),:,1:size(dp{i},3))=indexdp{i};
    n=n+size(dp{i},1);
end



function [d,xd,yd,indexd]=spheresort_current(lat,lon,lato,lono,cutoff)
dlat_cutoff=jrad2deg(cutoff./radearth);

sd=cell(length(lato),length(lono));
indexd=cell(length(lato),length(lono));
d=cell(length(lato),length(lono));
xd=cell(length(lato),length(lono));
yd=cell(length(lato),length(lono));

for j=1:length(lato)
     disp(['SPHERESORT computing latitude band ' int2str(j) ' of ' int2str(length(lato)) '.'])
     index=find(abs(lat-lato(j))<=dlat_cutoff);
     if ~isempty(index)
                  
         [latj,lonj]=vindex(lat,lon,index,1);
        
         lonomat=vrep(lono,length(index),1);
         latomat=lato(j)+zeros(size(lonomat));
                  
         [latjmat,lonjmat,indexmat]=vrep(latj,lonj,index,length(lono),2); 
         maxlat=min(abs(lato(j))+frac(360,2*pi)*frac(cutoff,radearth),90);
        
         %Form a ``wedge'' of nearby points; speeds things up quite a bit
         %dlon=abs(angle(exp(sqrt(-1)*frac(2*pi,360)*(lonjmat-lonomat))));
         
         %This is the same, but about twice as fast
         dlon=mod(lonjmat-lonomat,360);
         dlon=frac(2*pi,360)*min(dlon,360-dlon);
    
         %This comes from estimating distance only using delta-lon, and
         %re-arranging haversine formula with delta-lat set to zero in such
         %a way that trig functions are not applied to matrices
   
         fact=min(1,frac(sin(frac(cutoff,2*radearth)).^2,cosd(maxlat).*cosd(lato(j))));
         chi=2*asin(sqrt(fact));
         %if isnan(chi),chi=inf;end
         neari=find(dlon(:)<=chi);
           
         %length(find(neari))
         %dlon
         %Formerly I used this wedge, but it ends up being broader
         %This gives me the largest latitude where I might find
         %a point inside the search radius, but never > 90 degrees
         %neari=(dlon.*radearth.*cosd(maxlat))<=cutoff;
         %length(find(neari))
         
         %I had also tried using haversine to find exact distance, 
         %rearranging trig functions, but this ended up being slower
         
         [dj,xjmat,yjmat]=vzeros(size(latomat),'nan');
         [xjmat(neari),yjmat(neari),dj(neari)]=latlon2xy(latjmat(neari),lonjmat(neari),latomat(neari),lonomat(neari));
         dj(dj>cutoff)=nan;
         
         %Sort and return 
         if size(dj,1)>1
             [dj,sorter]=sort(dj);
             %After sorting the columns, I then rearrange all data
             %in each column according to the sort for that column
             
             sorterjj=vrep(1:size(sorter,2),size(sorter,1),1);
             indexsorter=sub2ind(size(xjmat),sorter,sorterjj);
             
             xjmat=xjmat(indexsorter);
             yjmat=yjmat(indexsorter);
             indexmat=indexmat(indexsorter);
         end
     
         bool=(sum(~isnan(dj),2)~=0);%Find all rows that have at least one non-nan
         vindex(dj,lonjmat,latjmat,indexmat,bool,1);
         
         for i=1:length(lono)
             ii=isfinite(dj(:,i));
             if ~isempty(ii)
                 d{j,i}=dj(ii,i);
                 xd{j,i}=xjmat(ii,i);
                 yd{j,i}=yjmat(ii,i);
                 indexd{j,i}=indexmat(ii,i);
             end
         end
     end
end

function [d,xd,yd,indexd]=spheresort_former(lat,lon,lato,lono,cutoff)

dlat_cutoff=jrad2deg(cutoff./radearth);

sd=cell(length(lato),length(lono));
indexd=cell(length(lato),length(lono));
d=cell(length(lato),length(lono));
xd=cell(length(lato),length(lono));
yd=cell(length(lato),length(lono));

[x,y,z]=latlon2xyz(lat,lon);

for j=1:length(lato)
     disp(['SPHERESORT computing latitude band ' int2str(j) ' of ' int2str(length(lato)) '.'])
     index=find(abs(lat-lato(j))<=dlat_cutoff);
     if ~isempty(index)
                  
         [xj,yj,zj,latj,lonj]=vindex(x,y,z,lat,lon,index,1);
         [xo,yo,zo]=latlon2xyz(lato(j),lono);
         %vsize(xo,yo,zo)
         [xomat,yomat,zomat,lonomat]=vrep(xo,yo,zo,lono,length(index),1);
         latomat=lato(j)+0*lonomat;
         
         [xjmat,yjmat,zjmat,latjmat,lonjmat,indexmat]=vrep(xj,yj,zj,latj,lonj,index,length(lono),2); 
     
         %vsize(xomat,xjmat,yomat,yjmat,zomat,zjmat)
         djtilde=sqrt(abs(squared(xomat-xjmat)+squared(yomat-yjmat)+squared(zomat-zjmat)));
         dj=2*radearth*asin(frac(1,2)*djtilde./radearth);  %asin gives phi in radians

         %cutofftilde=2*radearth*sin(frac(1,4*pi)*frac(cutoff,radearth));%=2*R*sin(phi/2)
         %neari=(djtilde<=cutofftilde);
         neari=(dj<=cutoff);
         [xtjmat,ytjmat]=vzeros(size(latomat),'nan');
         
         if anyany(neari)  %If any are true from neari
             [lattemp,lontemp]=xyz2latlon(xjmat(neari),yjmat(neari),zjmat(neari));
             %size(lattemp)
             %size(neari)
             %[xtjmat(neari),ytjmat(neari)]=latlon2xy(latomat(neari),lonomat(neari),lattemp,lontemp);
             
             [xtjmat(neari),ytjmat(neari)]=latlon2xy(lattemp,lontemp,latomat(neari),lonomat(neari));
             %[xtjmat(neari),ytjmat(neari)]=xyz2hor(latomat(neari),lonomat(neari),xjmat(neari),yjmat(neari),zjmat(neari));
         end

         dj(dj>cutoff)=nan;
         
         %Sort and return 
         if size(dj,1)>1
             [dj,sorter]=sort(dj);
             %After sorting the columns, I then rearrange all data
             %in each column according to the sort for that column
             
             sorterjj=vrep(1:size(sorter,2),size(sorter,1),1);
             indexsorter=sub2ind(size(xtjmat),sorter,sorterjj);
             
             xtjmat=xtjmat(indexsorter);
             ytjmat=ytjmat(indexsorter);
             indexmat=indexmat(indexsorter);
         end
       
         last=find(~isnan(vsum(dj,2)),1,'last');%Find last nan by summing over columns
         if ~isempty(last)
             vindex(dj,lonjmat,latjmat,indexmat,1:last,1);

             for i=1:length(lono)
                   ii=find(isfinite(dj(:,i)));
                   %length(ii)
                   if ~isempty(ii)
                       d{j,i}=dj(ii,i);
                       xd{j,i}=xtjmat(ii,i);
                       yd{j,i}=ytjmat(ii,i);
                       indexd{j,i}=indexmat(ii,i);
                   end
             end
         end
     end
end


function[]=spheresort_test


N=100000;
lat=rand(N,1)*180-90;
lon=rand(N,1)*360;

lono=(0:5:360);
lato=(-80:5:80);

cutoff=1000;

tic;[d,xd,yd,latd,lond]=spheresort(lat,lon,lat,lon,lato,lono,cutoff);etime1=toc;
%tic;[do,xdo,ydo,latdo,londo]=spheresort(lat,lon,lat,lon,lato,lono,cutoff);etime2=toc;

tic;
d2=zeros(size(d));
for i=1:length(lato)
    for j=1:length(lono)
        d2(i,j,:)=spheredist(latd(i,j,:),lond(i,j,:),lato(i),lono(j));
    end
end
etime2=toc;

%disp(['SPHERESORT was ' num2str(etime2./etime1) ' times faster than a simple loop.'])

bool1=zeros(length(lato),length(lono));
bool2=zeros(length(lato),length(lono));
for i=1:length(lato)
    for j=1:length(lono)
        bool1(i,j)=all(d(i,j,:)<=cutoff|isnan(d(i,j,:)));
        bool2(i,j)=aresame(d(i,j,:),d2(i,j,:),1e-6);
    end
end

reporttest('SPHERESORT all distances less than or equal to cutoff',allall(bool1))
reporttest('SPHERESORT verify distances',allall(bool2))

tic;[do,xdo,ydo,latdo,londo]=spheresort(lat,lon,lat,lon,lato,lono,cutoff,'former');etime2=toc;

tol=1e-8;
bool=[aresame(d,do,tol) aresame(xd,xdo,tol) aresame(yd,ydo,tol) aresame(latd,latdo,tol) aresame(lond,londo,tol)];
disp(['SPHERESORT current algorithm was ' num2str(etime2./etime1) ' times faster than former algorithm.'])
reporttest('SPHERESORT two algorithm versions match',allall(bool))

% ii=15;jj=30;
% plot(squeeze(lond(ii,jj,:)),squeeze(latd(ii,jj,:)),'.')
% hold on
% plot(lono(jj),lato(ii),'r*') 
% [lat,lon]=xy2latlon(squeeze(xd(ii,jj,:)),squeeze(yd(ii,jj,:)),lato(ii),lono(jj));
% plot(lon,lat,'go')


 
function[]=spheresort_parallel_test

%N=1000000;
N=10000;
lat=rand(N,1)*180-90;
lon=rand(N,1)*360;
 
lono=(0:5:360);
lato=(-80:5:80);

cutoff=1000;

%tic;[do,xdo,ydo,latd,lond]=spheresort(lat,lon,lat,lon,lato,lono,cutoff);etime1=toc;
tic;[do,xdo,ydo,latdo,londo]=spheresort(lat,lon,lat,lon,lato,lono,cutoff);etime1=toc;
tic;[d,xd,yd,latd,lond]=spheresort(lat,lon,lat,lon,lato,lono,cutoff,'parallel',12);etime2=toc;
bool=aresame(d,do)&&aresame(xd,xdo)&&aresame(yd,ydo)&&aresame(latd,latdo)&&aresame(lond,londo);
reporttest('SPHERESORT standard and parallel algorithms match',allall(bool))

%Wtih 10,000 datapoints parallel is 1/2 as fast
%With 100,000 datapoints parallel is 2x as fast 
%With 1,000,000 data points is 4x as fast
disp(['SPHERESORT parallel algorithm was ' num2str(etime1./etime2) ' times faster than standard algorithm.'])


%\************************************************************************


lato=[-90:0.1:90]';
cutoff=[100 200 400 800];
[chi,fact,maxlat]=vzeros(length(lato),length(cutoff));
for i=1:length(cutoff)
    for j=1:length(lato)
        maxlat(j,i)=min(abs(lato(j))+frac(360,2*pi)*frac(cutoff(i),radearth),90);
        fact(j,i)=min(1,frac(sin(frac(cutoff(i),2*radearth)).^2,cosd(maxlat(j,i)).*cosd(lato(j))));
        chi(j,i)=2*asin(sqrt(fact(j,i)));
    end
end
figure,plot(chi*360/2/pi,lato),ylim([-90 90])
