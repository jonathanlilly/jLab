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
%   The field VARNAME must have three dimensions, latitude, longitude, and 
%   time.  By default, the dimension order is expected to be (1) latitude,
%   (2) longitude, and (3) time.  NCINTERP(...,'transpose') specifies that
%   the order is instead (1) longitude, (2) latitude, and (3) time.
%
%   The NetCDF fil must also have variables named 'lat', 'lon', and 'num'
%   which give the values along the dimensions.  Both lat and lon are 
%   oriented in order of increasing values (i.e. lon cannot be 
%   discontinuous), while the time NUM is given in Matlab's DATENUM format.   
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
%   (C) 2016--2021 J.M. Lilly --- type 'help jlab_license' for details
 
filename=varargin{1};
numi=varargin{2};
lati=varargin{3};
loni=varargin{4};

if ~ischar(varargin{end})
    M=varargin{end};
    varargin=varargin(1:end-1);
else
    M=16;
end

str='fas';
orderstr='reg';
for i=1:2
    if ischar(varargin{end})
        if length(varargin{end})>=3
            if strcmp(varargin{end}(1:3),'loo')||strcmp(varargin{end}(1:3),'fas')
                str=varargin{end};
                varargin=varargin(1:end-1);
            elseif strcmp(varargin{end}(1:3),'reg')||strcmp(varargin{end}(1:3),'tra')
                orderstr=varargin{end};
                varargin=varargin(1:end-1);
            end
        end
    end
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

%determine whether to use +/- 180 or 360
if anyany(lon<0)
    loni=deg180(loni);
elseif anyany(lon>180)
    loni=deg360(loni);
end

bperiodic=false;
if (max(lon)-min(lon)+(lon(2)-lon(1)))==360
    disp('NCINTERP detecting 360 degrees of longitude.')
    %[min(lon) max(lon)]
    lon=[lon(end)-360;lon;lon(1)+360];
    %[min(lon) max(lon)]
    bperiodic=true;
end
num=double(ncread(filename,'num'));


%datestr(minmin(numi)),datestr(maxmax(numi))
%datestr(minmin(num)),datestr(maxmax(num))

%Set dates to NaNs when outside of range
bool=(numi<=minmin(num))|(numi>=maxmax(num));
numi(bool)=nan;

%size(numi)
%minmin(numi),maxmax(numi)
%length(find(~isnan(numi(:))))

if length(find(~isnan(numi(:))))>1
    a=find(num>=minmin(numi),1,'first');  %Note that a>=2
    b=find(num<=maxmax(numi),1,'last')+1; %Note that b<=length(num) 
else
    a=[];
    b=[];
end

if isempty(a)||isempty(b)
    for i=1:length(varnames)
        varargout{i}=nan*numi;
    end
else
    if strcmpi(orderstr(1:3),'reg')
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
    elseif  strcmpi(orderstr(1:3),'tra')
        for i=1:length(varnames)
            disp(['NCINTERP interpolating variable ' int2str(i) ' of ' int2str(length(varnames)) '.'])
            if strcmp(str(1:3),'loo')
                temp=ncinterp_one_loop_lonlat(filename,num,lat,lon,numi,lati,loni,varnames{i},blatup,blatdown,bperiodic,a,b);
            elseif strcmp(str(1:3),'fas')
                temp=ncinterp_one_fast_lonlat(filename,num,lat,lon,numi,lati,loni,varnames{i},blatup,blatdown,bperiodic,a,b,M);
            end
            if bwascell
                lattemp=lati;
                col2cell(lattemp,temp);
            end
            varargout{i}=temp;
        end
    end
end

function[xk]=ncinterp_one_fast_lonlat(filename,num,lat,lon,numi,lati,loni,varname,blatup,blatdown,bperiodic,a,b,M)
maxN=ceil(M*frac(1000*1000*1000,length(lat)*length(lon)*8)); %maximum number of time slices
xk=nan*zeros(size(numi));
n=0;

%a,N,b,return

for i=a:maxN:b
    n=n+1;
    disp(['NCINTERP reading chunk ' int2str(n) ' of ' int2str(length(a:maxN:b)) '.'])
    if n==1
        %Initialize xp
        xp=double(ncread(filename,varname,[1 1 1],[inf inf 1]));
        if blatdown,xp=[xp(:,1) xp];end    %Extra column for southern latitudes
        if blatup,  xp=[xp xp(:,end)];end  %Extra column for northern latitudes
    else
        if bperiodic
            xp=x([2:end-1],:,end);%Set xp back to last page of x
        else
            xp=x(:,:,end);
        end
    end
    
    %a,b,maxN,length(num)-i+1,b-i+1
    pagecount=min([maxN,length(num)-i+1,b-i+1]);
    %pagecount
    
    %datestr(num(i))
    %datestr(minmin(numi))
    %datestr(num(i+pagecount-1))
    %datestr(maxmax(numi))
    
    xn=double(ncread(filename,varname,[1 1 i],[inf inf pagecount]));
    
    x=vzeros(length(lon),length(lat),size(xn,3)+1);
    jj=1+blatdown:length(lat)-blatup;
    
    if bperiodic
        x([1 end],jj,1)=xp([end 1],jj);      %Extra row for longitudes
        x([1 end],jj,2:end)=xn([end 1],:,:); %Extra row for longitudes
        x(2:end-1,jj,1)=xp(:,jj);            %Central portion
        x(2:end-1,jj,2:end)=xn;              %Central portion
    else
        x(:,jj,1)=xp(:,jj);
        x(:,jj,2:end)=xn;
    end
    
    if blatdown,x(:,1,:)=x(:,2,:);end         %Extra column for southern latitudes
    if blatup,  x(:,end,:)=x(:,end-1,:);end   %Extra column for northern latitudes
    
    %vsize(x,xn,xp)
    
    %vsize(num,lat,lon,numi,lati,loni,x,xn,xp)%return
    bool=(numi>=num(i))&(numi<=num(i+pagecount-1));
    %vsize(numi,lati,loni,bool)
    [numk,latk,lonk]=vindex(numi(:),lati(:),loni(:),bool(:),1);  
    %round([min(lon) min(lonk) max(lon) max(lonk)])
    %vsize(lat,lon,num(i-1:i+pagecount-1),x,latk,lonk,numk)
    xk(bool) = interp3(lat,lon,num(i-1:i+pagecount-1),x,latk,lonk,numk,'linear');
end

%length(find(isfinite(xk)))

bool=isnan(xk)&~isnan(lati);
xk(bool)=inf;

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
    
    pagecount=min([maxN,length(num)-i+1,b-i+1]);
    xn=double(ncread(filename,varname,[1 1 i],[inf inf pagecount]));
    
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
    
    bool=(numi>=num(i))&(numi<=num(i+pagecount-1));
    [numk,latk,lonk]=vindex(numi,lati,loni,bool,1);  
%    vsize(lon,lat,num((i-1):(i+Ni)),x,lonk,latk,numk)
    %Note that INTERP3 bizarrely has reversed the first two input arguments
%    minmin(diff(lon)),minmin(diff(lat)),minmin(diff(num((i-1):(i+Ni))))
%     length(find(bool))
%     [1,length(find(isfinite(x)))]
%     plot(lon)
    xk(bool) = interp3(lon,lat,num(i-1:i+pagecount-1),x,lonk,latk,numk,'linear');
end

bool=isnan(xk)&~isnan(lati);
xk(bool)=inf;

function[xk]=ncinterp_one_loop_lonlat(filename,num,lat,lon,numi,lati,loni,varname,blatup,blatdown,bperiodic,a,b)

xn=double(ncread(filename,varname,[1 1 1],[inf inf 1]));
xk=nan*zeros(size(numi));
for i=a:b
    disp(['NCINTERP reading time ' int2str(i-a+1) ' of ' int2str(b-a+1) '.'])
    xp=xn;
    xn=double(ncread(filename,varname,[1 1 i],[inf inf 1]));
    
    x=vzeros(length(lon),length(lat),2);
    %ii=1+blatdown:size(xp,1)+blatdown;
    jj=1+blatdown:length(lat)-blatup;
    
    if bperiodic   
        x([1 end],jj,1)=xp([end 1],:); %Extra rows for longitudes
        x([1 end],jj,2)=xn([end 1],:); %Extra rows for longitudes
        x(2:end-1,jj,1)=xp;            %Central portion
        x(2:end-1,jj,2)=xn;            %Central portion
    else
        x(:,jj,1)=xp;
        x(:,jj,2)=xn;
    end
    
    if blatdown
        x(:,1,:)=x(:,2,:);         %Extra column for southern latitudes
    end
    if blatup
        x(:,end,:)=x(:,end-1,:);   %Extra column for northern latitudes
    end
    
    bool=(numi>=num(i-1))&(numi<=num(i));
    [numk,latk,lonk]=vindex(numi,lati,loni,bool,1);  
    %Note that INTERP3 bizarrely has reversed the first two input arguments
    xk(bool) = interp3(lat,lon,[num(i-1) num(i)],x,latk,lonk,numk,'linear');
end

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



