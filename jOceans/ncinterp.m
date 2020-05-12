function[varargout]=ncinterp(varargin)
%NCINTERP  One-line interpolation from 3D lat/lon/time field in NetCDF file.
%
%   X=NCINTERP(FILENAME,NUM,LAT,LON,VARNAME) interpolates the field VARNAME
%   from the NetCDF file FILENAME onto positions given by date numbers NUM,
%   latitudes LAT, and longitudes LON, and returns the results in X.
%
%   The interpolation is linear interpolation performed with INTERP3.
%   Longitude boundaries and near-polar latitudes are both accommodated.
%
%   The NetCDF file must a have a particular format.  The field VARNAME has
%   three dimensions, latitude, longitude, and time, in that order.  The 
%   file also has variables named 'lat', 'lon', and 'num' which give the 
%   values along the dimensions.  Both lat and lon are oriented in order of
%   increasing values.  The time NUM is given in Matlab's DATENUM format.   
%
%   NCINTERP detects if longitude in FILENAME spans 360 degrees, and if so,
%   interpolation across the longitudinal boundaries is correctly handled.
%   Note that longitudes can be specified as [-180,180] or [0,360].
%
%   NUM, LAT, and LON are all arrays of the same size, or else cell arrays
%   of numeric arrays.  In the latter case, X will also be a cell array.
%
%   [X1,X2,...,XN]=NCINTERP(FILENAME,NUM,LAT,LON,NAME1,NAME2,...,NAMEN),
%   with multiple variable names input, also works, for example:
%
%        [uwnd,vwnd]=ncinterp(filename,num,lat,lon,'uwnd','vwnd');
%
%   FILENAME should specify a file on the Matlab search path, or should
%   include the full pathname, e.g. FILENAME='/Volumes/Data/ncep_winds.nc'.  
%
%   Failures in the interpolation, for example, if the data falls outside 
%   the region covered by the NetCDF file, are set to values of INF.
%   __________________________________________________________________
%
%   Algorithms
%  
%   NCINTERP has two different approaches to reading the data.  The default
%   behavior is to load in time chunks such that each double-precision 
%   variable read in is smaller than a specified size, chosen as 16 GB.
%
%   NCINTERP(...,M) instead uses M gigabytes as the cutoff for the chunk. 
%   This method tends to be relatively fast.
%
%   NCINTERP(...,'loop') alternatively reads each time slice invididually
%   This method is very slow and is recommended only for testing purposes.
%   __________________________________________________________________
%
%   Handling domain edges
%
%   If the NetCDF file spans 360 degrees of longitude, the first and last
%   longitude are periodically wrapped to provide continuous interpolation.
%
%   To prevent the interpolation from failing very near the poles, the
%   first and last latitudes are also repeated, provided the last latitude
%   is within one latitude spacing interval of the pole. 
%   __________________________________________________________________
%
%   Usage: x=ncinterp(filename,num,lat,lon,varname);
%          [x1,x2]=ncinterp(filename,num,lat,lon,varname1,varname2);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2016--2018 J.M. Lilly --- type 'help jlab_license' for details
 
filename=varargin{1};
numi=varargin{2};
lati=varargin{3};
loni=varargin{4};

str='fas';

if ~ischar(varargin{end})
    M=varargin{end};
    varargin=varargin(1:end-1);
else
    M=16;
end

if ischar(varargin{end})
    if length(varargin{end})>=3
        if strcmp(varargin{end}(1:3),'loo')||strcmp(varargin{end}(1:3),'fas')
            str=varargin{end};
            varargin=varargin(1:end-1);
        end
    end
else
    str='loop';
end


varnames=varargin(5:end);

if iscell(numi)
    bwascell=true;
    [numi,lati,loni]=cell2col(numi,lati,loni);
else
    bwascell=false;
end


lat=double(ncread(filename,'lat'));

dlat=lat(end)-lat(end-1);
blatup=(lat(end)~=90)&&((lat(end)+dlat>=90));
if blatup
    lat=[lat;90];
end

dlat=lat(2)-lat(1);
blatdown=(lat(1)~=-90)&&((lat(1)-dlat<=-90));
if blatdown
    lat=[-90;lat];
end

lon=double(ncread(filename,'lon'));

%blatup,blatdown
%lonold=lon;
%lon=deg360(ncread(filename,'lon'));
%if anyany(abs(diff(lon))>180)
%    lon=lonold;
%end
    
bperiodic=false;
if (max(lon)-min(lon)+(lon(2)-lon(1)))==360
    disp('NCINTERP detecting 360 degrees of longitude.')
    lon=[lon(end)-360;lon;lon(1)+360];
    bperiodic=true;
end
num=double(ncread(filename,'num'));

%determine whether to use +/- 180 or 360
if anyany(lon<0)
    loni=deg180(loni);
elseif anyany(lon>180)
    loni=deg360(loni);
end

%Set dates to NaNs when outside of range
bool=(numi<=min(num))|(numi>=max(num));
numi(bool)=nan;

a=find(num>=min(numi),1,'first');  %Note that a>=2
b=find(num<=max(numi),1,'last')+1; %Note that b<=length(num) 

for i=1:length(varnames)
    disp(['NCINTERP interpolating variable ' int2str(i) ' of ' int2str(length(varnames)) '.'])
    if strcmp(str(1:3),'loo')
        temp=ncinterp_one_loop(filename,num,lat,lon,numi,lati,loni,varnames{i},blatup,blatdown,bperiodic,a,b);
    elseif strcmp(str(1:3),'fas')
        temp=ncinterp_one_fast(filename,num,lat,lon,numi,lati,loni,varnames{i},blatup,blatdown,bperiodic,a,b,M);
    end
    if bwascell
        lattemp=lati;
        col2cell(lattemp,temp);
    end
    varargout{i}=temp;
end

function[xk]=ncinterp_one_fast(filename,num,lat,lon,numi,lati,loni,varname,blatup,blatdown,bperiodic,a,b,M)
N=ceil(M*frac(1000*1000*1000,length(lat)*length(lon)*8));
   
% size(lon)
% lon(1:10)
% minmin(diff(lon))
% figure,plot(lon)
% return

xk=nan*zeros(size(numi));
n=0;
for i=a:N:b
    n=n+1;
    disp(['NCINTERP reading chunk ' int2str(n) ' of ' int2str(length(a:N:b)) '.'])
    if n==1
        %Initialize xp
        xp=double(ncread(filename,varname,[1 1 1],[inf inf 1]));
        if blatdown,xp=[xp(1,:);xp];end    %Extra row for southern latitudes
        if blatup,  xp=[xp;xp(end,:)];end  %Extra row for northern latitudes
    else
        if bperiodic
            xp=x(:,[2:end-1],end);%Set xp back to last page of x
        else
            xp=x(:,:,end);
        end
    end
    xn=double(ncread(filename,varname,[1 1 i],[inf inf min(N,length(num)-i+1)]));
    
    x=vzeros(length(lat),length(lon),size(xn,3)+1);
    ii=1+blatdown:length(lat)-blatup;
    
    %vsize(num,lat,lon,numi,lati,loni,x,xn,xp),N
    if bperiodic
        x(ii,[1 end],1)=xp(ii,[end 1]);      %Extra columns for longitudes
        x(ii,[1 end],2:end)=xn(:,[end 1],:); %Extra columns for longitudes
        x(ii,2:end-1,1)=xp(ii,:);            %Central portion
        x(ii,2:end-1,2:end)=xn;              %Central portion
    else
        x(ii,:,1)=xp(ii,:);
        x(ii,:,2:end)=xn;
    end
    
    
    if blatdown,x(1,:,:)=x(2,:,:);end         %Extra row for southern latitudes
    if blatup,  x(end,:,:)=x(end-1,:,:);end   %Extra row for northern latitudes
    
    Ni=min(N-1,length(num)-i);
    bool=(numi>=num(i-1))&(numi<=num(i+Ni));
    [numk,latk,lonk]=vindex(numi,lati,loni,bool,1);  
%    vsize(lon,lat,num((i-1):(i+Ni)),x,lonk,latk,numk)
    %Note that INTERP3 bizarrely has reversed the first two input arguments
%    minmin(diff(lon)),minmin(diff(lat)),minmin(diff(num((i-1):(i+Ni))))
%     length(find(bool))
%     [1,length(find(isfinite(x)))]
%     plot(lon)
    xk(bool) = interp3(lon,lat,num((i-1):(i+Ni)),x,lonk,latk,numk,'linear');
end

%length(find(isfinite(xk)))

bool=isnan(xk)&~isnan(lati);
xk(bool)=inf;

function[xk]=ncinterp_one_loop(filename,num,lat,lon,numi,lati,loni,varname,blatup,blatdown,bperiodic,a,b)

xn=double(ncread(filename,varname,[1 1 1],[inf inf 1]));
xk=nan*zeros(size(numi));
for i=a:b
    disp(['NCINTERP reading time ' int2str(i-a+1) ' of ' int2str(b-a+1) '.'])
    xp=xn;
    xn=double(ncread(filename,varname,[1 1 i],[inf inf 1]));
    
    x=vzeros(length(lat),length(lon),2);
    %ii=1+blatdown:size(xp,1)+blatdown;
    ii=1+blatdown:length(lat)-blatup;
    
    if bperiodic   
        x(ii,[1 end],1)=xp(:,[end 1]); %Extra columns for longitudes
        x(ii,[1 end],2)=xn(:,[end 1]); %Extra columns for longitudes
        x(ii,2:end-1,1)=xp;            %Central portion
        x(ii,2:end-1,2)=xn;            %Central portion
    else
        x(ii,:,1)=xp;
        x(ii,:,2)=xn;
    end
    
    if blatdown
        x(1,:,:)=x(2,:,:);         %Extra row for southern latitudes
    end
    if blatup
        x(end,:,:)=x(end-1,:,:);   %Extra row for northern latitudes
    end
    
    bool=(numi>=num(i-1))&(numi<=num(i));
    [numk,latk,lonk]=vindex(numi,lati,loni,bool,1);  
    %Note that INTERP3 bizarrely has reversed the first two input arguments
    xk(bool) = interp3(lon,lat,[num(i-1) num(i)],x,lonk,latk,numk,'linear');
end

bool=isnan(xk)&~isnan(lati);
xk(bool)=inf;

%tic;uwnd1=ncinterp('gom_ccmp_winds.nc',num,lat,lon,'uwnd');toc
%tic;uwnd2=ncinterp('gom_ccmp_winds.nc',num,lat,lon,'uwnd','loop');toc
%aresame(uwnd1,uwnd2,1e-8) 
%Fast version takes 2 seconds vs. 3 minutes for loop



