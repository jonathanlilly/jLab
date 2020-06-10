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
%   XS and YS are the corresponding coordinates of the sorted data points 
%   in a local tangent plane about each grid point. 
% 
%   LAT and LON are arrays of the same size into data point locations.
%   Any NaN values in LAT/LON are ignored.  
%
%   LATO and LONO are arrays of length M and N, say, specifying the
%   latitudes and longitudes of an M x N matrix of grid points, i.e.
%
%       LATO =  [LATO_1;    LONO= [LONO_1 LONO_2 ... LONO_N]. 
%                LATO_2; 
%                  ...
%                LATO_M]
%
%   SPHERESORT only computes distances for nearby points, up to a distance 
%   of CUTOFF, which is given in kilometers.  CUTOFF should be less than 
%   1/4 of the circumference of the earth, RADEARTH * PI/2, so that sorted
%   (LAT,LON) points will lie in the same hemisphere as the grid points.
%
%   The output arrays are M numerical arrays arranged as a length M cell
%   array.  That is, there is one cell per latitude band. Each numerical 
%   array has N columns, i.e., the number of longitudes, with the number of
%   rows varying between arrays.  
%
%   DS gives the distances of data points less than CUTOFF kilometers 
%   from the (m,n)th grid point, sorted in order of increasing distance.  
%   Entries farther than CUTOFF in all output fields are filled with NaNs.
%
%   XS and YS give the coordinates, in kilometers, of the corresponding 
%   data points in a plane tangent to the (m,n)th grid point. 
%
%   The choice to put latitude bands into cell arrays is made because the
%   number of points less than CUTOFF distance tends to vary with latitude,
%   and for convenience in parallelizing POLYSMOOTH.
%   _________________________________________________________________
% 
%   Limiting output dimension
%
%   SPHERESORT(LAT,LON,LATO,LONO,[CUTOFF JMAX]), where the fifth input 
%   argument is a 2-vector, additionally specifies that number of rows of
%   in each cell of the output will be no larger than JMAX.  This option is
%   useful for the 'fixed population' algorithm in POLYSMOOTH.
%   _________________________________________________________________
% 
%   Masking out non-mapping points
%
%   SPHERESORT(...,'mask',BOOL,...), where BOOL is a boolean array with
%   LENGTH(LATO) rows and LENGTH(LONO), only returns DS, XS, and YS at
%   points for which BOOL is true.  Where BOOL is false, the output arrays
%   will consist of all NaNs.
%
%   This is useful in preventing SPHERESORT from sorting data for locations
%   one has no intention of mapping.  For example, in creating a map of an
%   oceanic quantity, BOOL would be set to false for all land grid points.
%   _________________________________________________________________
% 
%   Additional input parameters
%
%   Let's say some additional variables Z1, Z2,...,ZK are given at
%   the data locations LAT, LON.  Then 
%   
%       [DS,XS,YS,Z1S,Z2S,...,ZKS]=
%           SPHERESORT(LAT,LON,Z1,Z2,...,ZK,LATO,LONO,CUTOFF);
%
%   also returns the sorted values of these variables.
%
%   Z1D, Z2D,...,ZKD are the same size as the other output arguments, and 
%   give the values of Z1, Z2,...,ZK sorted according to distance.
%   _________________________________________________________________
% 
%   Outputting an index
%
%   [DS,XS,YS,INDEXS]=SPHERESORT(LAT,LON,LATO,LONO,CUTOFF), with no
%   additional input arguments, returns a cell array of indices, INDEXS,
%   which is the same size as the other output fields.  This can be used 
%   together with POLYSMOOTH_EXTRACT to find sorted field values.  
%
%   If there are many fields to be mapped all using the same grid, this
%   approach is preferred to directly inputting Z1, Z2,...,ZK as in the 
%   previous step.  Instead the ZK would be output by POLYSMOOTH_EXTRACT.
%
%   For details, see the "One grid, many fields" section in POLYSMOOTH.  
%   _________________________________________________________________
%  
%   Parellel algorithm
%
%   SPHERESORT(...'parallel') parallelizes the computation using a parfor
%   loop over latitudes, which can speed things up dramatically.  This 
%   requires Matlab's Parallel Computing Toolbox to be installed.
%
%   SPHERESORT will then using an existing parallel pool, or if one does 
%   not exist, a pool will be created using all availabale workers.
%
%   SPHERESORT(...'parallel',Nworkers) alternately specifies the number of
%   workers to use. If you run into memory constraints, reduce Nworkers.
%   _________________________________________________________________
%
%   See also TWODSORT, POLYSMOOTH.
%
%   'spheresort --t' runs a test.
%
%   Usage: [ds,xs,ys,indexs]=spheresort(lat,lon,lato,lono,cutoff);
%          [ds,xs,ys]=spheresort(lat,lon,lato,lono,cutoff,'mask',bool);
%          [ds,xs,ys,z1s,z2s]=spheresort(lat,lon,z1,z2,lato,lono,cutoff);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2020 J.M. Lilly --- type 'help jlab_license' for details
 

if strcmpi(varargin{1}, '--t')
    spheresort_parallel_test;
    return
end

lat=varargin{1};
lon=varargin{2};

str='serial';
bool=[];
Nworkers=[];

for i=1:2
    if ischar(varargin{end-1})&&~ischar(varargin{end})
        if strcmpi(varargin{end-1}(1:3),'mas')
            bool=varargin{end};
        elseif strcmpi(varargin{end-1}(1:3),'par')
            str=varargin{end-1};
            Nworkers=varargin{end};
        end
        varargin=varargin(1:end-2);
    else
        if ischar(varargin{end})
             str=varargin{end};
             varargin=varargin(1:end-1);
        end
    end
end


pool = gcp('nocreate');
if strcmpi(str(1:3),'par')
    if isempty(Nworkers)
        if isempty(pool)
            parpool('local');
        end
    else
        if ~isempty(pool)
            parpool('local',Nworkers);
        elseif pool.NumWorkders~=Nworkers
            parpool('local',Nworkers);
        end
    end
end
    
lato=varargin{end-2};
lono=varargin{end-1};
cutoff=varargin{end};
args=varargin(3:end-3);

%In case max # points is not input
Ncutoff=inf;
if length(cutoff)==2
    Ncutoff=cutoff(2);
    cutoff=cutoff(1);
end

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

%Don't both sorting nan values
indexo=find(isfinite(lat)&isfinite(lon));

if ~isempty(indexo)
    vcolon(lat,lon);
    vindex(lat,lon,indexo,1);
else
    disp(['No finite data values.']), return
end

vcolon(lato,lono);
lono=lono';

if isempty(bool)
    bool=true(length(lato),length(lono));
end
 

if strcmpi(str(1:3),'par')
    [d,xd,yd,indexd]=spheresort_parallel(lat,lon,lato,lono,bool,cutoff,Ncutoff);
else
    [d,xd,yd,indexd]=spheresort_current(lat,lon,lato,lono,bool,cutoff,Ncutoff,[]);
end

varargout{1}=d;
varargout{2}=xd;
varargout{3}=yd;

for j=1:length(indexd)
    nonnani=~isnan(indexd{j});
    indexd{j}(nonnani)=indexo(indexd{j}(nonnani));
end

%Finally, account for extra input arguments
% for k=1:length(args)
%     temp=d;
%     for j=1:length(d)
%         nonnani=~isnan(indexd{j});
%         temp{j}(nonnani)=args{k}(indexd{j}(nonnani));
%     end
%     varargout{k+3}=temp;
% end

for k=1:length(args)
    varargout{k+3}=polysmooth_sortfield(indexd,args{k});
end

varargout{length(args)+4}=indexd;

disp('SPHERESORT finished.')



function [d,xd,yd,indexd]=spheresort_parallel(lat,lon,lato,lono,bool,cutoff,Ncutoff)

%Parallelize by looping over latitude
d=cell(length(lato),1);
xd=cell(length(lato),1);
yd=cell(length(lato),1);
indexd=cell(length(lato),1);

parfor i=1:length(lato)
    [d1,xd1,yd1,indexd1]=spheresort_current(lat,lon,lato(i),lono,bool(i,:),cutoff,Ncutoff,[i length(lato)]);
    d{i}=d1{1};
    xd{i}=xd1{1};
    yd{i}=yd1{1};
    indexd{i}=indexd1{1};
end

%delete(gcp)

function [d,xd,yd,indexd]=spheresort_current(lat,lon,lato,lono,bool,cutoff,Ncutoff,J)
dlat_cutoff=jrad2deg(cutoff./radearth);

d=cell(length(lato),1);
xd=cell(length(lato),1);
yd=cell(length(lato),1);
indexd=cell(length(lato),1);
for j=1:length(lato)
     if isempty(J)
         disp(['SPHERESORT computing latitude band ' int2str(j) ' of ' int2str(length(lato)) '.'])
     else
         disp(['SPHERESORT computing latitude band ' int2str(J(1)) ' of ' int2str(J(2)) '.'])
     end
     
     %first we remove land (or non-mapping) points     
     oceanindex=find(bool(j,:));
     lono1=lono(oceanindex);

     index=find(abs(lat-lato(j))<=dlat_cutoff);
     if ~isempty(index)
                  
         [latj,lonj]=vindex(lat,lon,index,1);
        
         lonomat=vrep(lono1,length(index),1);
         latomat=lato(j)+zeros(size(lonomat));
                  
         [latjmat,lonjmat,indexmat]=vrep(latj,lonj,index,length(lono1),2); 
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
         xjmat(isnan(dj))=nan;
         yjmat(isnan(dj))=nan;
         indexmat(isnan(dj))=nan;

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
     
         bool1=(sum(~isnan(dj),2)~=0);%Find all rows that have at least one non-nan
         vindex(dj,xjmat,yjmat,indexmat,bool1,1);
         if size(dj,1)>Ncutoff
             vindex(dj,xjmat,yjmat,indexmat,1:Ncutoff,1); 
         end
         %there was previously a major bug here where I did not index xjmat
         %and yjmat in the previous call
         
         d{j}(1:size(dj,1),1:length(lono))=nan;
         xd{j}(1:size(dj,1),1:length(lono))=nan;
         yd{j}(1:size(dj,1),1:length(lono))=nan;
         indexd{j}(1:size(dj,1),1:length(lono))=nan;
         
         d{j}(:,oceanindex)=dj;
         xd{j}(:,oceanindex)=xjmat;
         yd{j}(:,oceanindex)=yjmat;
         indexd{j}(:,oceanindex)=indexmat;
     end
end


function[]=spheresort_parallel_test

%N=1000000;
N=1000000;
lat=rand(N,1)*180-90;
lon=rand(N,1)*360;
 
lono=(0:5:360);
lato=(-80:5:80);

cutoff=1000;

%tic;[do,xdo,ydo,latd,lond]=spheresort(lat,lon,lat,lon,lato,lono,cutoff);etime1=toc;
tic;[d,xd,yd,latd,lond]=spheresort(lat,lon,lat,lon,lato,lono,cutoff,'parallel');etime2=toc;
tic;[do,xdo,ydo,latdo,londo]=spheresort(lat,lon,lat,lon,lato,lono,cutoff);etime1=toc;
bool=aresame(d,do)&&aresame(xd,xdo)&&aresame(yd,ydo)&&aresame(latd,latdo)&&aresame(lond,londo);
reporttest('SPHERESORT standard and parallel algorithms match',allall(bool))

%older spmd algorithm:
%With 10,000 datapoints parallel is 1/2 as fast
%With 100,000 datapoints parallel is 2x as fast 
%With 1,000,000 data points is 4x as fast

%new parfor algorithm:
%With 10,000 datapoints parallel is 1.9x as fast
%With 100,000 datapoints parallel is 2.1x as fast 
%With 1,000,000 data points is 2.5x as fast

disp(['SPHERESORT parallel algorithm was ' num2str(etime1./etime2) ' times faster than standard algorithm.'])

%%to plot the nth slice of d
% ds1=nan*ssh;
% for i=1:length(ds)
%      if ~isempty(ds{i})
%          for j=1:length(lon)
%              ds1(i,j)=ds{i}(1,j);
%          end
%      end
% end


