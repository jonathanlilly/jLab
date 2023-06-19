function[znew]=interplatlon(varargin)
%INTERPLATLON  Interpolation for working with latitude and longitude.
%
%   INTERPLATLON performs 2D or 3D interpolation appropriate for latitude
%   and longitude.  It deals with issues arising from the periodicity of 
%   longitude in addition to performing several other convenient functions. 
%   __________________________________________________________________
%
%   2D Interpolation
%
%   Z=INTERPLATLON(LATO,LONO,[],ZO,LAT,LON) performs a 2D interpolation
%   from the array ZO, which has latitude varying in LENGTH(LATO) rows and
%   longitude varying in LENGTH(LONO) columns, onto locations (LAT,LON).
%
%   Here LAT and LON are 2D arrays of the same size.  LATO and LONO should
%   1D arrays such that ZO is plaid; otherwise see SPHEREINTERP.
%
%   A difficulty in interpolating latitude and longitude is that longitude
%   is periodic, for example, if LONO=0.5:1:359.5.  Then the interpolation
%   can fail if 0<LON<0.5 or 359.5<LON<360.  INTERPLATLON checks for LONO 
%   being periodic, and if so it handles the interpolation appropriately by
%   appending the last column of ZO to its front and the last to is back.
%
%   INTERPLATLON accounts for the fact that LONO may have an arbitrary 
%   starting value, not necessarily -180 or 0, as its initial point, e.g.
%   LONO may begin at -20 and end at 340.  Note that both LONO and LON are
%   interpreted as angles on the complex plane, so the modifications
%   LONO+n*360 or LON+n*360 for integer n would not change the result.  
%
%   ZO may have multiple elements along its 3rd, 4th, or 5th dimension.  Z
%   will then have its first two dimensions matching LAT and LON, and
%   higher dimensions matching ZO.  
%
%   Z=INTERPLATLON(LATO,LONO,[],ZO,LAT,LON,[],BOOL), where BOOL is a
%   boolean array with size SIZE(LAT), sets the corresponding entries of 
%   Z to NaNs.  BOOL could for example be true over land.  If Z has more
%   than 2 dimensions, this is done for each SIZE(LAT) element of Z.
%   __________________________________________________________________
%
%   Interpolation algorithm
%
%   By default, INTERPLATLON uses Matlab's INTERP2 with linear 
%   interpolation.  Alternatively, INTERPLATLON(...,METHOD) uses the method
%   specified by the string METHOD, as described in INTERP2.
% 
%   Data points for which the chosen interpolation method does not return a
%   valid value are then filled using nearest-neighbor interpolation.
%
%   INTERPLATON(...,'fill'), if LAT and LON are matrices of the from output 
%   by MESHGRID, will instead fill missing data points from the linear
%   interpolation with the results of a bin-averaging.  
%
%   This typically affects only a small number of points, but is useful in
%   controlling exactly what constitutes missing data.  
%
%   For example, when interpolating an ocean field from a finer grid to a 
%   coaser grid, linear interpolation will result in somewhat degraded 
%   continental boundaries. The 'fill' algorithm puts back in a few data 
%   points that would otherwise be missed by linear interpolation.
%   __________________________________________________________________
%
%   3D Interpolation
%
%   Z=INTERPLATLON(LATO,LONO,TO,ZO,LAT,LON,T) similarly performs a 3D
%   interpolation, where TO is now time, with SIZE(ZO,3) being the same as
%   LENGTH(TO).  LAT, LON, and T are all the same size.  
%
%   Again, INTERPLATLON(...,METHOD) specifies the interpolation method,
%   which defaults to 'linear'.
%
%   Note that the 'fill' option does not work with 3D interpolation.
%   __________________________________________________________________
%
%   Usage: z=interplatlon(lato,lono,[],zo,lat,lon);
%          z=interplatlon(lato,lono,[],zo,lat,lon,[],bool);
%          z=interplatlon(lato,lono,to,zo,lat,lon,t);
%          z=interplatlon(lato,lono,to,zo,lat,lon,t,bool);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018--2022 J.M. Lilly --- type 'help jlab_license' for details
 
str='nofill';
interpstr='linear';
for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'fil')||strcmpi(varargin{end}(1:3),'nof')
            str=varargin{end};
        else
            interpstr=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end
tol=1e-6;

lato=varargin{1}(:);
lono=varargin{2}(:);
to=varargin{3}(:);
zo=varargin{4};
lat=varargin{5};
lon=varargin{6};

if length(varargin)>=7
    t=varargin{7};
else 
    t=[];
end

if length(varargin)==8
    land=varargin{8};
else
    land=[];
end

dlono=lono(2)-lono(1);

if aresame(rot(2*pi*(lono(1)-dlono)/360),rot(2*pi*lono(end)/360),tol)
    if ~aresame(lono(end)-lono(end-1),dlono,tol)
        error('INTERPLATLON detecting periodicity, but LONO spacing is not the same at the beginning and end.')
    end 
    disp('INTERPLATLON detects periodic longitude array; wrapping ends.')
    lono=[lono(1)-dlono;lono;lono(end)+dlono];
    %this appends the last column to the front and the first to the end
    zo=vshift(zo,-1,2);
    zo(:,end+1,:,:,:)=zo(:,1,:,:,:);
    zo(:,end+1,:,:,:)=zo(:,2,:,:,:);    
end

%[min(lono) max(lono)]

%lono(1)
%account for the fact that lono may not start at -180 or 0
%This sets the minimum value of the longitude I'm interpolating onto to be
%halfway between the first two longitudes we're interpolating from
lon=frac(360,2*pi)*angle(rot(frac(2*pi,360)*(lon-lono(1)-dlono/2-180)))+lono(1)+dlono/2+180;
%xx=frac(360,2*pi)*angle(rot(frac(2*pi,360)*([-180:.1:180]-25-1/2-180)))+25+180+1/2;
%min(xx)
%minmin(lon),maxmax(lon),minmin(lono),maxmax(lono)
%minmin(lon)>minmin(lono)&&maxmax(lon)<maxmax(lono)

ismeshgrid=false;
if strcmpi(str(1:3),'fil')
    %is what we're interpolating onto a valid meshgrid?
    b1=~isscalar(lon);
    b2=aresame(vrep(lon(1,:),size(lon,1),1),lon,tol);
    b3=aresame(vrep(lat(:,1),size(lat,2),2),lat,tol);
    b4=allall(abs(diff(lon)-(lon(2)-lon(1)))<tol);
    b5=allall(abs(diff(lat)-(lat(2)-lat(1)))<tol);
    ismeshgrid=b1&&b2&&b3&&b4&&b5;
    %[b1,b2,b3,b4,b5]
    if ismeshgrid
        disp('INTERPLATLON detects MESHGRID format output; filling gaps with bin averages.')
        dlat=lat(2,1)-lat(1,1);
        dlon=lon(1,2)-lon(1,1);
        %lon(1,1),lon(1,2),dlon,dlat
        xedges=lon(1)-dlon/2:dlon:lon(end)+dlon/2;
        yedges=lat(1)-dlat/2:dlat:lat(end)+dlat/2;
    else
        disp('INTERPLATLON not detecting MESHGRID format, so not filling gaps with bin averages.')
    end
end
%[long,latg]=meshgrid(lon,lat);      %this is grid we're interpolating onto


if isempty(t)
    %2D interpolation
    znew=nan*zeros([size(lat,1) size(lat,2) size(zo,3) size(zo,4) size(zo,5)]);
    for i=1:size(zo,3)
        for j=1:size(zo,4)
            for k=1:size(zo,5)
                %vsize(lonog,latog,zo(:,:,i,j,k),lon,lat)
                temp=interp2(lono,lato,zo(:,:,i,j,k),lon,lat,interpstr);
                if ismeshgrid
                    %vsize(lonog,latog,zo(:,:,i,j,k),lon,lat,xedges,yedges)
                    %fill missing data from linear interpolation with bin averages
                    [lonog,latog]=meshgrid(lono,lato);  %this is the grid we interpolating from
                    [temp2,lon2,lat2]=twodstats(lonog,latog,zo(:,:,i,j,k),xedges,yedges);
                    temp(isnan(temp))=temp2(isnan(temp));
                    %aresame(vcolon(lon(1,:)),lon2(:),tol)
                    %aresame(vcolon(lat(:,1)),lat2(:),tol)
                else
                    %fill missing data from linear interpolation with nearest neighbor
                    bool=isnan(temp);
                    temp(bool)=interp2(lono,lato,zo(:,:,i,j,k),lon(bool),lat(bool),'nearest');
                end
                
                if ~isempty(land)
                    temp(land)=nan;
                end
                temp(isnan(lat))=nan;
                znew(:,:,i,j,k)=temp;
            end
        end
    end
else
    %3D interpolation
    znew=nan*zeros([size(lat,1) size(lat,2) size(zo,4) size(zo,5)]);
    for j=1:size(zo,4)
        for k=1:size(zo,5)
            %vsize(lonog,latog,zo(:,:,i,j,k),lon,lat)
            temp=interp3(lono,lato,to,zo(:,:,:,j,k),lon,lat,t,interpstr);
            bool=isnan(temp);
            temp(bool)=interp3(lono,lato,to,zo(:,:,:,j,k),lon(bool),lat(bool),t(bool),'nearest');
            if ~isempty(land)
                temp(land)=nan;
            end
            temp(isnan(lat))=nan;
            znew(:,:,j,k)=temp;
        end
    end
end

