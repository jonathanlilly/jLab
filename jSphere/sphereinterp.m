function[varargout]=sphereinterp(varargin)
%SPHEREINTERP  Fast linear interpolation on the sphere from non-plaid grids.
%   _______________________________________________________________________
%   
%   *|* sphereinterp1.png --- Result of mapping with SPHEREINTERP.  With 
%   JLAB installed, type 'jhelp sphereinterp' to view this image. *|*
%   _______________________________________________________________________
%
%   SPHEREINTERP takes a 2D array ZO, defined for latitudes and longitudes
%   LATO and LONO, and linearly interpolates the values of ZO onto a plaid 
%   grid defined by 1-D arrays LAT and LON, leading to the output array Z.
%
%   LATO, LONO, and ZO are all the same size.  It is not necessary for LATO
%   and LONO to be plaid, i.e. created by MESHGRID. LATO should increase 
%   monotonically along array *rows*, while LONO should increase
%   monotonically along array *columns*.  Both should vary smoothly as a
%   function of the underlying array coordinates.
%
%   SPHEREINTERP is called twice with different arguments, first with four 
%   arguments and then with five, as follows:   
%
%      [DX,DY,INDEX,BOOL]=SPHEREINTERP(LATO,LONO,LAT,LON);
%      Z=SPHEREINTERP(DX,DY,INDEX,BOOL,ZO);
%
%   The first phase computes indices into the original LATO, LONO arrays 
%   for the target LAT, LON values.  The second phase uses these indices to 
%   find linearly interpolated values through an array lookup.  
%
%   The first phase is computationally expensive, but the second is very
%   fast.  For intepolating many different fields on the same grid, one 
%   only needs to make the first call to SPHEREINTERP one time, thus
%   giving a computationally efficient way to perform the interpolation.
%
%   If LATO and LONO are a plaid grid, INTERP2 will be much faster and
%   should be used instead.
%
%   It is not necessary to understand the details of these two calls to
%   SPHEREINTERP, however for completeness they are described below. 
%   _________________________________________________________________
%
%   Periodicity in longitude
%
%   SPHEREINTERP(LATO,LONO,LAT,LON,'periodic') adapts the first phase 
%   calculation to account for periodicity in longitude.  This option 
%   should be used if LATO and LONO represent the entire sphere. 
%
%   In this case, LATO and LONO should be periodic in the column direction,
%   that is, the column just to the right of the last column is understood
%   to be represented by the first column. 
%   ___________________________________________________________________
%
%   Parallelization
%
%   SPHEREINTERP(LATO,LONO,LAT,LON,'parallel') parallelizes the first 
%   phase of the calculation using a PARFOR loop. This requires that 
%   Matlab's Parallel Computing Toolbox be installed. 
%   ___________________________________________________________________
%
%   Condition number
%
%   [DX,DY,INDEX,BOOL,C]=SPHEREINTERP(LATO,LONO,LAT,LON) with five output
%   arguments for the first phase call also outputs C, the condition number
%   of the Jacobian matrix of the LATO,LONO grid at each LAT/LON location.  
%
%   If there are locations where this matrix is ill-conditioned, one may 
%   form a threshold on C to revert to the nearest-neighbor interpolation, 
%   which can be output as described below under 'Second phase details.'
%   __________________________________________________________________
%
%   First phase details: index computation
%
%   [DX,DY,INDEX,BOOL]=SPHEREINTERP(LATO,LONO,LAT,LON) finds the point
%   in the original LATO, LONO fields nearest each target LAT, LON point, 
%   as well as the three points adjacent to this closest point.  
%   
%   DX and DY are arrays of the same size as LAT and LON giving the column
%   and row deviations, respectively, within LATO and LONO from the closest
%   point in those arrays to the linearly interpolated LAT/LON value.
%
%   The closest original point to each target point is found using 
%   SPHEREDIST.  The DX and DY arrays are found from the observed lat/lon
%   deviations by inverting (using MATINV) the Jacobian matrix describing
%   the variation of LATO and LONO with respect to the array coordinates.
%
%   INDEX is a cell array with four elements, each of which is the same
%   size as LAT and LON.  It gives the locations within LATO and LONO of
%   the closest point to each target point, and of three adjacent points.
%
%   Specifically, INDEX{1} gives the index into LATO and LONO of the
%   original point nearest each target point.  INDEX{2} is the index into
%   the closer of the two original points in an adjacent *column* to the
%   closest point, INDEX{3} is likewise in an adjacent *row* to the closest 
%   point, and INDEX{4} is in both an adjacent row and an adjacent column.
%
%   BOOL is a boolean array that is true whenever all four of the indices 
%   in INDEX are well-defined. 
%   _______________________________________________________________________
%
%   Second phase details: linear interpolation
%
%   Z=SPHEREINTERP(DX,DY,INDEX,BOOL,ZO) then computes Z as a weighted sum 
%   of the ZO values from the four points identified in the first stage:
%
%       Z1 = (1-ABS(DX)).* (1-ABS(DY)) .* ZO(INDEX{1});
%       Z2 =    ABS(DX) .* (1-ABS(DY)) .* ZO(INDEX{2});
%       Z3 = (1-ABS(DX)).* ABS(DY)     .* ZO(INDEX{3});
%       Z4 =    ABS(DX) .* ABS(DY)     .* ZO(INDEX{4});
%       Z  = Z1 + Z2 + Z3 + Z4;
%
%   For the case in which LATO and LONO are a plaid grid, this matches the
%   results from Matlab's INTERP2 to numerical precision. 
%
%   At any points where the linear interpolation fails or yields NaNs---as 
%   can happen where the inversion of the Jacobian is not numerically 
%   stable, or where interpolated locations point to undefined values of 
%   the original field (e.g. near coastlines)---the nearest-neighbor fit
%   is substituted in place of the linear interpolation.
%
%   [Z,Z1]=SPHEREINTERP(DX,DY,INDEX,BOOL,ZO) also outputs Z1, the nearest-
%   neighbor interpolation, for comparison.  Then LENGTH(FIND(Z==Z1)) gives
%   the number of points for which the nearest-neighbor fit has been used.
%   _______________________________________________________________________
%   
%   *|* sphereinterp2.png ---Example of using SPHEREINTERP.  With 
%   JLAB installed, type 'jhelp sphereinterp' to view this image. *|*
%   _______________________________________________________________________
%
%   Example
%
%   The above two figures provide an example of using SPHEREINTERP.  The
%   relevant code is contained in MAKEFIGS_SPHEREINTERP.
%
%   The left-hand figure shows sea surface height from an ocean model which
%   uses a "tripolar" grid.  Note that at high Northern latitudes, features 
%   are not only strectched, they are also distorted.  Using SPHEREINTERP 
%   the model fields are mapped onto a regular lat/lon grid, seen at right. 
%
%   The figure at the top of the page is the sea surface height gradient 
%   magnitude of the interpolated field, showing a very high level of
%   detail.  Some minor artifacts at two locations in the Arctic reveal
%   the locations where the model grid is singular.  
%
%   Computing the initial mapping coefficients on the models rougly 1/8 
%   degree grid is computationally very expensive, and takes about 1.5
%   hours on a 12 core Mac Pro working in parallel mode.  After that,
%   however, the mapping for each time step takes only about one second
%   using the second call to SPHEREINTERP.  
%
%   The model data is from a simulation using the GFDL's Generalized Ocean
%   Layered Model, or GOLD, kindly provided by Harper Simmons at the
%   University of Alaska Fairbanks.
%   _________________________________________________________________
%
%   See also SPHERESORT, JSPHERE.
%
%   'sphereinterp --t' runs a test.
%   'sphereinterp --f' generates the two figures shown above.
%
%   Usage: [dx,dy,index,bool]=sphereinterp(lato,lono,lat,lon,'parallel');
%          z=sphereinterp(dx,dy,index,bool,zo);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017--2022 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    sphereinterp_test,return
elseif strcmp(varargin{1}, '--f')
    makefigs_sphereinterp,return
end

wrapstr='endpoint';
parstr='serial';
for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
            parstr=varargin{end};
            varargin=varargin(1:end-1);
%         elseif strcmpi(varargin{end}(1:3),'pre')||strcmpi(varargin{end}(1:3),'int')
%             str=varargin{end};
%             varargin=varargin(1:end-1);
         elseif strcmpi(varargin{end}(1:3),'per')||strcmpi(varargin{end}(1:3),'end')
            wrapstr=varargin{end};
            varargin=varargin(1:end-1);
        end
    end
end


if length(varargin)==4
    lato=varargin{1};
    lono=varargin{2};
    lat=varargin{3};
    lon=varargin{4};
    [dx,dy,index,bool,C]=sphereinterp_prep(lato,lono,lat,lon,parstr,wrapstr);
    
    varargout{1}=dx;
    varargout{2}=dy;
    varargout{3}=index;
    varargout{4}=bool;
    varargout{5}=C;
else
    dx=varargin{1};
    dy=varargin{2};
    index=varargin{3};
    bool=varargin{4};
    zo=varargin{5};
    
    [z,z1,z2,z3,z4]=vzeros(size(dx),'nan');
    z(bool)=zo(index{1}(bool));
    
    z1(bool)=(1-abs(dx(bool))).*(1-abs(dy(bool))).*zo(index{1}(bool));
    z2(bool)=abs(dx(bool)).*(1-abs(dy(bool))).*zo(index{2}(bool));
    z3(bool)=(1-abs(dx(bool))).*abs(dy(bool)).*zo(index{3}(bool));
    z4(bool)=abs(dx(bool)).*abs(dy(bool)).*zo(index{4}(bool));
    znew=z1+z2+z3+z4;
    
    %Revert to nearest if linear fails
    bool=~isfinite(znew)&isfinite(z);
    znew(bool)=z(bool);

    %Put nans where nearest-neighbor is a nan
    znew(isnan(z))=nan;
    varargout{1}=znew;
    varargout{2}=z;
end


function[dx,dy,index,bool,C]=sphereinterp_prep(lato,lono,lat,lon,parstr,wrapstr)

if ~aresame(size(lato),size(lono))
    [lono,lato]=meshgrid(lono,lato);
end

%--------------------------------------------------------------------------
%Preparing for linear interpolation by finding indices
if ~aresame(size(lat),size(lon))
    [long,latg]=meshgrid(lon,lat);
    [ii,jj]=vzeros(size(latg),'nan');
else
    [ii,jj]=vzeros(size(lat),'nan');
end

%This version works when I have plaid (meshgrid) output
[jjo,iio]=meshgrid(1:size(lato,2),1:size(lato,1));
if strcmpi(parstr(1:3),'par')
    parfor i=1:length(lat)
        disp(['SPHEREINTERP finding indices for latitude band ' int2str(i) ' of ' int2str(length(lat)) '.'])
        %Sort in latitude bands
        bool=(vshift(lato,-1,1)<=lat(i))&(vshift(lato,1,1)>lat(i));
        index=find(bool);
        [ii1,jj1]=vzeros(size(lon),'nan');
        for j=1:length(lon)
            d=spheredist(lato(index),lono(index),lat(i),lon(j));
            [mind,id]=min(d);
            %i,length(index),mind
            if ~isempty(id)
                ii1(j)=iio(index(id));
                jj1(j)=jjo(index(id));
            end
        end
        ii(i,:)=ii1;
        jj(i,:)=jj1;
    end
else
    %This is just an exact copy of the above, but without the parallel loop
    for i=1:length(lat)
        disp(['SPHEREINTERP finding indices for latitude band ' int2str(i) ' of ' int2str(length(lat)) '.'])
        %Sort in latitude bands
        bool=(vshift(lato,-1,1)<=lat(i))&(vshift(lato,1,1)>lat(i));
        index=find(bool);
        [ii1,jj1]=vzeros(size(lon),'nan');
        for j=1:length(lon)
            d=spheredist(lato(index),lono(index),lat(i),lon(j));
            [mind,id]=min(d);
            %i,length(index),mind
            if ~isempty(id)
                ii1(j)=iio(index(id));
                jj1(j)=jjo(index(id));
            end
        end
        ii(i,:)=ii1;
        jj(i,:)=jj1;
    end
end
%--------------------------------------------------------------------------

%Create an index that looks up the locations of (lat,lon) within (lato,lono)
%This is the same size as (latg,long)
index=nan*zeros(size(ii));
bool=isfinite(ii)&isfinite(jj);
index(bool)=sub2ind(size(lato),ii(bool),jj(bool));
bool=isfinite(index);%Only those value of the index that are finite

% figure,plot(jj')
% minmin(ii)
% maxmax(ii)
% minmin(jj)
% maxmax(jj)
%length(find(~bool(:)))

%Nearest-neighbor lookup of z within zo
%z(bool)=zo(index(bool));

%Derivatives of lono and lato
dlondx=vdiff(unwrap(lono,2),2,wrapstr);
dlatdx=vdiff(lato,2,wrapstr);
dlondy=vdiff(unwrap(lono,1),1);
dlatdy=vdiff(lato,1);

dlondx=deg180(dlondx*2)/2;%Very important!  For the first central difference
%employed here, jumps across 360 tend to lead to *two* points exceeding
%180.  This sets those back to nearly zero, as expected.

%Values of derivatives of (lato,lono) observed at (lat,lon) locations
[M11,M12,M21,M22]=vzeros(size(latg),'nan');
M11(bool)=dlondx(index(bool));
M12(bool)=dlondy(index(bool));
M21(bool)=dlatdx(index(bool));
M22(bool)=dlatdy(index(bool));

clear M
M=zeros(2,2,size(latg,1),size(latg,2));
M(1,1,:,:)=M11;
M(1,2,:,:)=M12;
M(2,1,:,:)=M21;
M(2,2,:,:)=M22;

%M=zeros(size(latg,1),size(latg,2),2,2);
%M(:,:,1,1)=M11;
%M(:,:,1,2)=M12;
%M(:,:,2,1)=M21;
%M(:,:,2,2)=M22;

invM=matinv(M);%Inverting this matrix
%squeeze(M(10,10,:,:))*squeeze(invM(10,10,:,:))  %Just checking

%Values of nearest (lato,lono) observed at (lat,lon) locations
[loni,lati]=vzeros(size(latg),'nan');
loni(bool)=lono(index(bool));
lati(bool)=lato(index(bool));

%Difference in lat and lon from nearest-neighbor interpolated to grid values
dlat=latg-lati;
dlon=angle(rot(frac(2*pi,360)*(long-loni)))*frac(360,2*pi);

%dx=invM(:,:,1,1).*dlon+invM(:,:,1,2).*dlat;
%dy=invM(:,:,2,1).*dlon+invM(:,:,2,2).*dlat;

dx=squeeze(invM(1,1,:,:)).*dlon+squeeze(invM(1,2,:,:)).*dlat;
dy=squeeze(invM(2,1,:,:)).*dlon+squeeze(invM(2,2,:,:)).*dlat;

%Generally, these will be less than one, because they represent shifts
%within the (lato,lono) grid.  But a small fraction, for numerical reasons,
%tend to exceed one and these should just be set to zero. Setting the
%shifts to zero means we will use the nearest-neighbor fit instead.
isoutside=abs(dx)>1|abs(dy)>1;
dx(isoutside)=0;
dy(isoutside)=0;

%ii = *longitude* index, jj=*latitude* index
[jj,ii]=ind2sub(size(lato),index);

iishift=sign(dx);
jjshift=sign(dy);

%--------------------------------------------------------------------------
%Account for points on the edges of the domain 
if strcmpi(wrapstr(1:3),'per')
    bool=(ii+iishift==size(lato,2)+1);
    iishift(bool)=1-size(lato,2);
    bool=(ii+iishift==0);
    iishift(bool)=size(lato,2)-1;
else
    bool=(ii+iishift==size(lato,2)+1)|(ii+iishift==0);
    iishift(bool)=0;
end

bool=(jj+jjshift==size(lato,2)+1)|(jj+jjshift==0);
jjshift(bool)=size(lato,1);
%--------------------------------------------------------------------------

% minmin(ii+iishift)
% maxmax(ii+iishift)
% minmin(jj+jjshift)
% maxmax(jj+jjshift)
% size(lato)

%These should be very close to zero for test on meshgrid
%dlon2=nan*dlon;dlat2=nan*dlat;
%dlon2=(1-abs(dx)).*lono1(ii)+abs(dx).*lono1(ii+iishift)-long;
%dlat2=(1-abs(dy)).*lato1(jj)+abs(dy).*lato1(jj+jjshift)-latg;


clear index
index{1}=sub2ind(size(lato),jj,ii);                 %Shifting neither
index{2}=sub2ind(size(lato),jj,ii+iishift);         %Shifting lon
index{3}=sub2ind(size(lato),jj+jjshift,ii);         %Shifting lat
index{4}=sub2ind(size(lato),jj+jjshift,ii+iishift); %Shifting both
bool=isfinite(index{1}.*index{2}.*index{3}.*index{4});

%Condition number loop
C=nan*ones(size(M(:,:,1,1)));
for i=1:size(M,3)
    for j=1:size(M,4)
        if allall(isfinite(M(:,:,i,j)))
            C(i,j)=cond(squeeze(M(:,:,i,j)));
        end
    end
end


function[]=sphereinterp_test
 
load slasnapshot
lono=slasnapshot.lon;
lato=slasnapshot.lat;
zo=slasnapshot.sla;

%lon=-64:1/10:-46;lat=26:1/10:42;
lon=-64:1/5:-46;lat=26:1/5:42;

tic;[dx,dy,index,bool]=sphereinterp(lato,lono,lat,lon);toc
tic;[zlin,z2]=sphereinterp(dx,dy,index,bool,zo);toc

[long,latg]=meshgrid(lon,lat);
tic;znearest=interp2(lono,lato,zo,long,latg,'nearest');toc
tic;zlinear=interp2(lono,lato,zo,long,latg,'linear');toc

reporttest('SPHEREINTERP compare with INTERP2 for plaid grid',aresame(zlin,zlinear,1e-10))

%lon=-64:1/10:-46;lat=26:1/10:42;
lono=slasnapshot.lon;
lato=slasnapshot.lat;
zo=slasnapshot.sla;

lon=172:1/5:190;lat=26:1/5:42;

[dx,dy,index,bool]=sphereinterp(lato,lono,lat,lon,'parallel','periodic');
tic;[zlin,z2]=sphereinterp(dx,dy,index,bool,zo);toc

[long,latg]=meshgrid(lon,lat);
[lono,zo]=lonshift(0,lono,zo);

tic;znearest=interp2(lono,lato,zo,long,latg,'nearest');toc
tic;zlinear=interp2(lono,lato,zo,long,latg,'linear');toc

reporttest('SPHEREINTERP compare with INTERP2 for plaid grid, encompassing dateline',aresame(zlin,zlinear,1e-10))
