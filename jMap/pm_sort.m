function[varargout]=pm_sort(varargin)
%PM_SORT  Sort data points according to distance from grid points in 2D.
%
%   Given data points distributed in two dimensions along with a grid, 
%   PM_SORT returns the data points sorted by increasing distance from each
%   grid point, up to some specified maximum.  This works both for 
%   Cartesian geometry as well as on the surface of Earth. 
%
%   PM_SORT is called internally by POLYMAP.  However, for large problems
%   it may be preferable to call it externally, as documented in POLYMAP.
%   _________________________________________________________________
%   
%   Sorting on the plane
%
%   [DS,XS,YS]=PM_SORT(X,Y,XO,YO,CUTOFF) returns sorted distances D between
%   data points at locations X,Y and grid points at XO,YO.
%
%   X and Y are arrays of the same size into data point locations. Any non-
%   finite values in X and Y are ignored.
% 
%   XO and YO are arrays of length M and L, say, specifying the bin center
%   locations of an L x M matrix of grid points, i.e.
%
%      XO= [XO_1 XO_2 ... XO_M]      YO =  [YO_1;    
%                                           YO_2; 
%                                            ...
%                                           YO_L]
%
%   CUTOFF is the maximum distance to be included in the output arrays.
%
%   The output arrays are L numerical arrays arranged as a length L cell
%   array.  That is, there is one cell per element of YO. Each numerical 
%   array has M columns, i.e., the number of elements of XO, with the 
%   number of rows varying between arrays.  This organization was chosen to
%   in order to optimize memory usage and facilitate paralellization.
%
%   DS gives the distances SQRT((X-XO)^2+(Y-YO)^2) of all data points less
%   than the CUTOFF distance from the (l,m)th grid point, sorted in order
%   of increasing distance.  Entries farther than CUTOFF in all output
%   fields are filled with NaNs.
%
%   XS and YS are corresponding deviations X-XO and Y-YO from the grid 
%   point location to each data point.
%   _________________________________________________________________
% 
%   Limiting output dimension
%
%   [DS,XS,YS]=PM_SORT(X,Y,XO,YO,[CUTOFF NMAX]), where the fifth input 
%   argument is a 2-vector, additionally specifies that number of rows of
%   in each cell of the output will be no larger than NMAX.  This option is
%   useful for the 'fixed population' algorithm in POLYMAP.
%   _________________________________________________________________
% 
%   Additional input parameters
%
%   Let's say some additional variables Z1, Z2,...,ZK are given at the data
%   locations X,Y.  Then 
%   
%   [DS,XS,YS,Z1S,Z2S,...,ZKS]=
%
%                PM_SORT(X,Y,Z1,Z2,...,ZK,XO,YO,CUTOFF);
%
%   also returns the values of these variables.
%
%   Z1S, Z2S,...,ZKS are the same size as the other output arguments, and 
%   give the values of Z1, Z2,...,ZK sorted according to distance.
%   _________________________________________________________________
% 
%   Outputting an index
%
%   [DS,XS,YS,INDEXS]=PM_SORT(X,Y,XO,YO,CUTOFF), with no additional input 
%   arguments, returns a cell array of indices, INDEXS, which is the same 
%   size as the other output fields.  This can be used together with 
%   PM_INDEX to find sorted field values.  
%
%   If there are many fields to be mapped all using the same grid, this
%   approach is preferred to directly inputing Z1, Z2, ..., ZK as in the 
%   described previously. Instead the ZK would be output on as as-needed 
%   basis by PM_INDEX.
%
%   For usage details, see the "One grid, many fields" section in POLYMAP. 
%   _________________________________________________________________
%   
%   Sorting on the surface of the Earth
%
%   PM_SORT can also handle data points distributed on the Earth's surface.
%
%   [DS,XS,YS]=PM_SORT(LON,LAT,LONO,LATO,CUTOFF,'sphere') returns the great
%   circle distances DS between data points at locations LON, LAT and
%   nearby grid points located at LONO, LATO, sorted in order of distance. 
%   Longitude LON and latitude LAT are given in degrees.
%
%   XS and YS are the corresponding x- and y-coordinates, in kilometers, of
%   the sorted data points in a local tangent plane about each grid point. 
% 
%   The maximum value of CUTOFF is RADEARTH * PI/2, or a quarter of the 
%   circumference of the Earth, so that sorted (LAT,LON) points will 
%   necessarily lie in the same hemisphere as the grid points.  Larger 
%   values of CUTOFF will be set to this value.
%
%   As with sorting on the plane, given LATO of length L and LONO of length 
%   M, the output fields will again be L numerical arrays arranged as 
%   length L cell arrays, with each numerical array having M columns. 
%
%   PM_SORT uses lat/lon truncation criteria to exclude points that cannot 
%   be nearer than CUTOFF, avoiding unnecessary computations involving the 
%   great circle distance and considerably speeding up the sorting.
%
%   Note that since the number of points less than CUTOFF distance tends to 
%   vary strongly with latitude, the use of 3D arrays for the sorted fields 
%   would lead to large arrays with many NaNs.  The choice to put latitude
%   bands into cell arrays avoids this and thus uses less system memory. 
%
%   Note that all of the above options also work with the 'sphere' flag.
%   _________________________________________________________________
%   
%   Masking out regions on the Earth
%
%   When the 'sphere' flag is input, PM_SORT(...,'mask',BOOL,...), with
%   BOOL being an L x M matrix, will returns output fields only at points 
%   for which BOOL is true.  Where BOOL is false, the output arrays will 
%   consist of all NaNs.
%
%   This is useful in preventing PM_SORT from sorting data for locations
%   one has no intention of mapping.  For example, in creating a map of an
%   oceanic quantity, BOOL would be set to false for all land grid points.
%   _________________________________________________________________
%  
%   Other options
%
%   PM_SORT supports a number of other options, as documented more 
%   completely in POLYMAP.
% 
%   PM_SORT(...,CUTOFF,...,'population',N) specifies a fixed
%   population fit with population N and maximum bandwidth CUTOFF.  N can 
%   be a scalar or an L x M matrix.
%
%   PM_SORT(...,'verbose') displays the row it is currently working on.
%
%   PM_SORT(...,'parallel') parallelizes the computation using PARFOR.
%
%   PM_SORT(...'parallel',Nworkers) parallelizes the computation with a
%   specified number of workers.
%   _________________________________________________________________
% 
%   See also POLYMAP and PM_INDEX.
%
%   'pm_sort --t' runs a test.
%
%   Usage: [ds,xs,ys,indexs]=pm_sort(x,y,xo,yo,cutoff);
%          [ds,xs,ys,zs,indexs]=pm_sort(x,y,z,xo,yo,cutoff);
%          [ds,xs,ys,zs,indexs]=pm_sort(lon,lat,z,lon,lat,cutoff,'sphere');
%          [ds,xs,ys,z1s,z2s,...,zNs]=pm_sort(x,y,z1,z2,...,zN,xo,yo,cutoff);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    pm_sort_consistency_test 
    pm_sort_parallel_test
    return
end

%--------------------------------------------------------------------------
%sort out string input arguments
bool=[];
Nworkers=[];
parstr='serial';
geostr='cartesian';
verstr='quiet';
Npop=nan;

for i=1:6
    if ischar(varargin{end-1})&&~ischar(varargin{end})
        if strcmpi(varargin{end-1}(1:3),'mas')
            bool=varargin{end};
        elseif strcmpi(varargin{end-1}(1:3),'par')
            parstr=varargin{end-1};
            Nworkers=varargin{end};
        elseif strcmpi(varargin{end-1}(1:3),'pop')
            Npop=varargin{end};
        end
        varargin=varargin(1:end-2);
    elseif ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'ser')||strcmpi(varargin{end}(1:3),'par')
            parstr=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'car')||strcmpi(varargin{end}(1:3),'sph')
            geostr=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'ver')||strcmpi(varargin{end}(1:3),'qui')
            verstr=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

%--------------------------------------------------------------------------
%sort out rest of input arguments
xdata=varargin{1};
ydata=varargin{2};
xo=varargin{end-2}(:)';
yo=varargin{end-1}(:);

index=(1:length(xdata(:)))';%index into original data locations

if isempty(bool)
    bool=true(length(yo),length(xo));
end
 
cutoff=max(varargin{end});
varargin=varargin(3:end-3);

if ~aresame(size(xdata),size(ydata))
    error('The first two input arguments must be the same size.')
end

%Remove non-finite values
index_finite=find(isfinite(xdata)&isfinite(ydata));
xdata=xdata(index_finite);
ydata=ydata(index_finite);
index=index(index_finite);
if isempty(xdata)
     disp(['No finite data values.']), return
end

if numel(Npop)==1,   Npop=Npop+zeros(length(yo),length(xo));     end
%--------------------------------------------------------------------------
%begin parallelization if requested
pool = gcp('nocreate');
if strcmpi(parstr(1:3),'par')
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

%--------------------------------------------------------------------------
%implement loops 

L=length(yo);
ds=cell(L,1);
xs=cell(L,1);
ys=cell(L,1);
indexs=cell(L,1);

if strcmpi(geostr(1:3),'car')
    if strcmpi(parstr(1:3),'ser')
        for i=1:L
            if strcmpi(verstr(1:3),'ver')
                disp(['PM_SORT working on row ' int2str(i) ' of ' int2str(L) '.'])
            end
            [ds{i},xs{i},ys{i},indexs{i}]=pm_sort_cartesian(xdata,ydata,xo,yo(i),bool(i,:),index,cutoff,Npop(i,:));
        end
    elseif strcmpi(parstr(1:3),'par')
        parfor i=1:L
            if strcmpi(verstr(1:3),'ver')
                disp(['PM_SORT working on row ' int2str(i) ' of ' int2str(L) '.'])
            end
            [ds{i},xs{i},ys{i},indexs{i}]=pm_sort_cartesian(xdata,ydata,xo,yo(i),bool(i,:),index,cutoff,Npop(i,:));
            %Matlab doesn't like me to assign to varargout in the parfor loop
        end
    end
elseif strcmpi(geostr(1:3),'sph')
    if strcmpi(parstr(1:3),'ser')
        for i=1:L
            if strcmpi(verstr(1:3),'ver')
                disp(['PM_SORT working on row ' int2str(i) ' of ' int2str(L) '.'])
            end
            %note ydata=latitude, xdata=longitude, yo=latitude, xo=longitude
            [ds{i},xs{i},ys{i},indexs{i}]=pm_sort_sphere(xdata,ydata,xo,yo(i),bool(i,:),index,cutoff,Npop(i,:));
        end
    elseif strcmpi(parstr(1:3),'par')
        parfor i=1:L
            if strcmpi(verstr(1:3),'ver')
                disp(['PM_SORT working on row ' int2str(i) ' of ' int2str(L) '.'])
            end
            %note ydata=latitude, xdata=longitude, yo=latitude, xo=longitude
            [ds{i},xs{i},ys{i},indexs{i}]=pm_sort_sphere(xdata,ydata,xo,yo(i),bool(i,:),index,cutoff,Npop(i,:));
            %Matlab doesn't like me to assign to varargout in the parfor loop
        end
    end
end

varargout{1}=ds;
varargout{2}=xs;
varargout{3}=ys;

for i=1:length(varargin)
    varargout{i+3}=pm_index(indexs,varargin{i});
end

varargout{length(varargin)+4}=indexs;

function[ds,xs,ys,indexs]= pm_sort_cartesian(xdata,ydata,xo,yo,bool,index,cutoff,Npop)
%vsize(xdata,ydata,xo,yo)

Ncutoff=min(inf,max(Npop));

M=length(xo);
N=length(xdata);

oceanindex=find(bool);
xo=xo(oceanindex);

xp=xdata-xo;
yp=vrep(ydata-yo,length(xo),2);
index=vrep(index,length(xo),2);

d=sqrt(xp.^2+yp.^2);

d(d>cutoff)=nan;
xp(isnan(d))=nan;
yp(isnan(d))=nan;
index(isnan(d))=nan;

if ~allall(~isfinite(d))
    [dsort,sorter]=sort(d,'ascend');
    jj=vrep(1:length(xo),N,1);
    sorter=sub2ind(size(d),sorter,jj);

    xp=xp(sorter);
    yp=yp(sorter);
    index=index(sorter);

    Nmax=min(find(sum(isfinite(dsort),2),1,'last'),Ncutoff);

    %locations we're instructed not to map are set to nan
    ds=nan*zeros(Nmax,M);
    xs=nan*zeros(Nmax,M);
    ys=nan*zeros(Nmax,M);
    indexs=nan*zeros(Nmax,M);
    
    %vsize(ds,d)
    ds(:,oceanindex)=dsort(1:Nmax,:);
    xs(:,oceanindex)=xp(1:Nmax,:);
    ys(:,oceanindex)=yp(1:Nmax,:);
    indexs(:,oceanindex)=index(1:Nmax,:);
else
    ds=zeros(0,M);
    xs=zeros(0,M);
    ys=zeros(0,M);
    indexs=zeros(0,M);
end

function [ds,xs,ys,indexs]=pm_sort_sphere(lon,lat,lono,lato,bool,index,cutoff,Npop)

Ncutoff=min(inf,max(Npop));
cutoff=min(cutoff,radearth*pi/2);

dlat_cutoff=jrad2deg(cutoff./radearth);

%first we remove land (or non-mapping) points
oceanindex=find(bool);
lono1=lono(oceanindex);

%size(bool)
nearlat=find(abs(lat-lato)<=dlat_cutoff);

%vsize(nearlat,lat,lon,lato,dlat_cutoff)

M=length(lono);
if ~isempty(nearlat)

    [lat,lon,index]=vindex(lat,lon,index,nearlat,1);

    lonomat=vrep(lono1,length(nearlat),1);
    latomat=lato+zeros(size(lonomat));

    %vsize(lono1,lato,lonomat,latomat)

    [latmat,lonmat,indexmat]=vrep(lat,lon,index,length(lono1),2);
    maxlat=min(abs(lato)+frac(360,2*pi)*frac(cutoff,radearth),90);

    %Form a ``wedge'' of nearby points; speeds things up quite a bit
    %dlon=abs(angle(exp(sqrt(-1)*frac(2*pi,360)*(lonjmat-lonomat))));

    %This is the same, but about twice as fast
    %vsize(lonmat,lonomat)
    dlon=mod(lonmat-lonomat,360);
    dlon=frac(2*pi,360)*min(dlon,360-dlon);

    %This comes from estimating distance only using delta-lon, and
    %re-arranging haversine formula with delta-lat set to zero in such
    %a way that trig functions are not applied to matrices

    fact=min(1,frac(sin(frac(cutoff,2*radearth)).^2,cosd(maxlat).*cosd(lato)));
    chi=2*asin(sqrt(fact));
    %if isnan(chi),chi=inf;end
    nearlon=find(dlon(:)<=chi);

    %Formerly I used this wedge, but it ends up being broader
    %This gives me the largest latitude where I might find
    %a point inside the search radius, but never > 90 degrees
    %neari=(dlon.*radearth.*cosd(maxlat))<=cutoff;
    %length(find(neari))

    %I had also tried using haversine to find exact distance,
    %rearranging trig functions, but this ended up being slower

    [d,xmat,ymat]=vzeros(size(latomat),'nan');
    [xmat(nearlon),ymat(nearlon),d(nearlon)]=...
        latlon2xy(latmat(nearlon),lonmat(nearlon),latomat(nearlon),lonomat(nearlon));
  
    d(d>cutoff)=nan;
    xmat(isnan(d))=nan;
    ymat(isnan(d))=nan;
    indexmat(isnan(d))=nan;

    %Sort and return
    if size(d,1)>1
        [d,sorter]=sort(d);
        %After sorting the columns, I then rearrange all data
        %in each column according to the sort for that column

        sorteri=vrep(1:size(sorter,2),size(sorter,1),1);
        indexsorter=sub2ind(size(xmat),sorter,sorteri);

        xmat=xmat(indexsorter);
        ymat=ymat(indexsorter);
        indexmat=indexmat(indexsorter);
    end

%    bool1=(sum(~isnan(d),2)~=0);%Find all rows that have at least one non-nan
    bool1=(sum(isfinite(d),2)~=0);%Find all rows that have at least one finite value
    vindex(d,xmat,ymat,indexmat,bool1,1);
    if size(d,1)>Ncutoff
        vindex(d,xmat,ymat,indexmat,1:Ncutoff,1);
    end

    Nmax=size(d,1);

    %locations we're instructed not to map are set to nan
    ds=nan*zeros(Nmax,M);
    xs=nan*zeros(Nmax,M);
    ys=nan*zeros(Nmax,M);
    indexs=nan*zeros(Nmax,M);
    
    %vsize(ds,d)
    ds(:,oceanindex)=d;
    xs(:,oceanindex)=xmat;
    ys(:,oceanindex)=ymat;
    indexs(:,oceanindex)=indexmat;
else
    ds=zeros(0,M);
    xs=zeros(0,M);
    ys=zeros(0,M);
    indexs=zeros(0,M);
end

vswap(ds,xs,ys,indexs,inf,nan);

function[]=pm_sort_consistency_test

rng(1);
[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:200);

[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.125:3);
yo=(-3:.125:3);
[xg,yg]=meshgrid(xo,yo);

[ds,xs,ys,xs2,ys2]=pm_sort(xdata,ydata,xdata,ydata,xo,yo,10);
for i=1:length(xs)
    bool(i)=aresame(xs{i}+vrep(xo,size(xs{i},1),1),xs2{i})&aresame(ys{i}+yo(i),ys2{i});
end
reporttest('PM_SORT Cartesian case consistency, no missing data',allall(bool));


%Insert some NANs
xdata(1:7:end)=nan;
[ds,xs,ys,xs2,ys2]=pm_sort(xdata,ydata,xdata,ydata,xo,yo,10);
for i=1:length(xs)
    bool(i)=aresame(xs{i}+vrep(xo,size(xs{i},1),1),xs2{i})&aresame(ys{i}+yo(i),ys2{i});
end
reporttest('PM_SORT Cartesian case consistency, with missing data',allall(bool));


function[]=pm_sort_parallel_test

N=10000;
lat=rand(N,1)*180-90;
lon=rand(N,1)*360;
 
lono=(0:5:360);
lato=(-80:5:80);

cutoff=1000;

%tic;[do,xdo,ydo]=pm_sort(lat,lon,lato,lono,cutoff,'sphere');etime1=toc;
%tic;[d,xd,yd]=pm_sort(lat,lon,lato,lono,cutoff,'sphere','parallel');etime2=toc;
tic;[do,xdo,ydo,latdo,londo]=pm_sort(lat,lon,lat,lon,lato,lono,cutoff,'sphere');etime1=toc;
tic;[d,xd,yd,latd,lond]=pm_sort(lat,lon,lat,lon,lato,lono,cutoff,'sphere','parallel');etime2=toc;
bool=aresame(d,do)&&aresame(xd,xdo)&&aresame(yd,ydo)&&aresame(latd,latdo)&&aresame(lond,londo);
reporttest('PM_SORT with sphere flag, standard and parallel algorithms match',allall(bool))

%older spmd algorithm:
%With 10,000 datapoints parallel is 1/2 as fast
%With 100,000 datapoints parallel is 2x as fast 
%With 1,000,000 data points is 4x as fast

%previous parfor algorithm:
%With 10,000 datapoints parallel is 1.9x as fast
%With 100,000 datapoints parallel is 2.1x as fast 
%With 1,000,000 data points is 2.5x as fast

%current parfor algorithm:
%With 10,000 datapoints parallel is 3.6x as fast
%With 100,000 datapoints parallel is 2.6x as fast 
%With 1,000,000 data points is 2.4x as fast

disp(['PM_SORT with sphere flag, parallel algorithm was ' num2str(etime1./etime2) ' times faster than standard algorithm.'])

