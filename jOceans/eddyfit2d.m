function[varargout]=eddyfit2d(varargin)
%EDDYFIT2D  Least squares fit of 2D velocity data to an eddy profile.
%
%   FIT=EDDYFIT2D(NUM,LAT,LON,U,V,LATO,LONO,D,INDEX,Z0,ZA,ZB) performs a
%   least squares fit of 2D velocity measurements to a Guassian eddy.
%  
%   The fit involves four parameters, the eddy center location as well as
%   the core radius R and peak velocity V.
%   __________________________________________________________________
%
%   Input
%
%   The velocity measurements U and V are 3D arrays taken at times NUM,
%   latitudes LAT, and longitudes LON, each of which is a 1D array. U and V
%   have LENGTH(LAT) rows, LENGTH(LON) columns, and LENGTH(NUM) pages along
%   the third dimension.  NUM is in Matlab's DATENUM format.
%   
%   A fit is performed to all data points less than D kilometers distant
%   from the point with coordinates LATO and LONO, as well as over sets of 
%   time slices specified by INDEX.  
%
%   INDEX can be one of the following:
%
%        INDEX = empty, []  -- Fit to velocity averaged over all pages
%        INDEX = array   -- Fit to velocity averaged over INDEX pages
%        INDEX = cell of N arrays --  N fits, averaged over INDEX{1}, etc.
%
%   The arrays Z0, ZA, and ZB provide the initial guesses, lower bounds,
%   and upper bounds for three of the four parameters of the fit, with 
%
%        Z0 = [L0 R0 V0],  ZA = [LA RA VA],   ZB = [LB RB VB].
%
%   Here L is the distance of the eddy center to the origin (LATO,LONO) in 
%   kilometers, R is the eddy core radius in kilometers, and V is the
%   signed azimuthal velocity at the core radius R in cm/s.
%
%   To have no upper bound or lower bound, use e.g. VB = inf or VB = -inf.
%
%   Only three rather than four parameters are initialized becaues the fit
%   is performed in cylindical coordinates, with no bounds being set on the
%   value of the azimuth angle. 
%   __________________________________________________________________
%
%   Output
%
%   The output FIT is a structure with the following fields:
%
%        FIT.NUM        Date in DATENUM format, as input 
%        FIT.LONE       Estimated eddy center longitudes  
%        FIT.R          Estimated eddy core radii in kilometers
%        FIT.V          Estimated signed peak azimuthal velocity in cm/s
%        FIT.ERR        Total root-mean-square velocity error 
%        FIT.ERRP       Root-mean-square error from azimuthal velocities 
%        FIT.ERRN       Root-mean-square error from normal velocities 
%        FIT.FLAG       The exit flag from the optimization
%        FIT.COUNT      Total number of good velocity points in the fit
%        FIT.OUTPUT     Output fields from the optimization
%        FIT.LAT        Latitudes of observations, as input
%        FIT.LON        Longitudes of observations, as input 
%        FIT.UM         3D array of observed averaged zonal velocities 
%        FIT.VM         3D array of observed averaged meridional velocities 
%        FIT.UE         3D array of zonal velocities implied by eddy 
%        FIT.VE         3D array of meridional velocities implied by eddy 
%        
%   All fields except OUTPUT and those following are N x 1 arrays, where N 
%   is the total number of fits. N = 1 unless INDEX is a cell array, in 
%   which case N = LENGTH(INDEX).  
%
%   COUNT keeps track of the totalnumber of good data points involved in 
%   the fit. OUTPUT is an N x 1 cell array of structures containing
%   information from the optimatization routine. 
%
%   UM, VM, UE, and VE are of size LENGTH(LAT) x LENGTH(LON) x N, with UM 
%   and VM containing the averaged velocites observed over each element of
%   INDEX, and UE and VE containing the velocities implied by the eddy fit. 
%
%   The errors are related by ERR.^2=ERRP.^2+ERRN.^2.  The error is related
%   to the observed mean (MU,MV) and predicted velocities (UE,VE) by 
%   ERR.^2 = MEAN((UM-UE).^2+(VM-VE).^2), where MEAN is over all finite 
%   values of velocity in the time slices specified by INDEX. 
%   __________________________________________________________________
%
%   Parallelization
%
%   EDDYFIT2D(...,'parallel') parallelizes the computation with a PARFOR
%   loop over the elements of INDEX, provides this is a cell array.  This 
%   option requires Matlab's Parallel Computing Toolbox be installed.
%   __________________________________________________________________
%
%
%   Algorithm details
%
%   EDDYFIT2D works by applying FMINSEARCHBND to solve the constrained
%   optimization problem.  FMINSEARCHBND, by is an open-source modification 
%   of Matlab's standard FMINSEARCH, and is distributed with JLAB. 
%
%   In order to deal with the cyclic nature of the azimuth angle, two fits
%   are performed, first with bounds -pi and pi and guess 0, and second
%   with bounds 0 and 2pi and guess pi; then the better of these fits is
%   chosen.  This avoids the possibility of edge effects impacting the fit. 
%   __________________________________________________________________
%
%   Usage: fit=eddyfit2d(num,lat,lon,u,v,lato,lono,D,index,z0,za,zb);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    eddyfit2d_test,return
end
 
%index=[];
str='Gaussian';
parstr='ser';

tol=1e-6;

%if iscell(varargin{1})
%varargin=varargin(2:end);
%end

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'par')|| strcmpi(varargin{end}(1:3),'ser')
            parstr=varargin{end};
        else
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end
if length(varargin)==13
   tol=varargin{13};
end


num=varargin{1};
lat=varargin{2};
lon=varargin{3};
u=varargin{4};
v=varargin{5};
lato=varargin{6};
lono=varargin{7};
D=varargin{8};
index=varargin{9};
z0=varargin{10};
za=varargin{11};
zb=varargin{12};

[xg,yg]=meshgrid(lon,lat);
[x,y,d]=latlon2xy(yg,xg,lato,lono);
cv=u+1i*v;

r=x+1i*y;
bool=d<D;

options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',tol,'TolX',tol);

%z0=xguess(:);
%za=za(:)';
%zb=zb(:)';

if isempty(index)
    index=1:length(num);
end

if ~iscell(index)
    %r=vrep(r,length(num),3);
    %bool=vrep(bool,length(num),3);
    
    [late,lone,R,V,um,vm,ue,ve,err,errp,errn,count,flag,output]=...
        eddyfit2d_one(r,cv(:,:,index),bool,z0,za,zb,lato,lono,options,str);
else
    
    %r=vrep(r,length(num),3);
    %bool=vrep(bool,length(num),3);
    
    [late,lone,R,V,err,errp,errn,flag,count]=vzeros(length(index),1,nan);
    [um,vm,ue,ve]=vzeros(size(r,1),size(r,2),length(index),nan);
    output=cell(length(index),1);
    
    if strcmpi(parstr(1:3),'ser')
        for i=1:length(index)
            disp(['EDDYFIT2D computing fit ' int2str(i) ' of ' int2str(length(index)) '.'])
            [late(i),lone(i),R(i),V(i),um(:,:,i),vm(:,:,i),ue(:,:,i),ve(:,:,i),err(i),errp(i),errn(i),count(i),flag(i),output{i}]=...
                eddyfit2d_one(r,cv(:,:,index{i}),bool,z0,za,zb,lato,lono,options,str);
        end
    else
        parfor i=1:length(index)
            disp(['EDDYFIT2D computing fit ' int2str(i) ' of ' int2str(length(index)) '.'])
            [late(i),lone(i),R(i),V(i),um(:,:,i),vm(:,:,i),ue(:,:,i),ve(:,:,i),err(i),errp(i),errn(i),count(i),flag(i),output{i}]=...
                eddyfit2d_one(r,cv(:,:,index{i}),bool,z0,za,zb,lato,lono,options,str);
        end
    end
end

%make struct late lone R V ue ve err count flag struct

%lat=late;
%lon=lone;
%u=ue;
%v=ve;

struct.num=num;
struct.late=late;
struct.lone=lone;
struct.R=R;
struct.V=V;
struct.err=err;
struct.errp=errp;
struct.errn=errn;
struct.count=count;
struct.flag=flag;
struct.output=output;
struct.lat=lat;
struct.lon=lon;
struct.um=um;
struct.vm=vm;
struct.ue=ue;
struct.ve=ve;
varargout{1}=struct;


%make struct lat lon R V u v err count flag output
%varargout{1}=late;
%varargout{2}=lone;
%varargout{3}=R;
%varargout{4}=V;
%varargout{1}=ue;
%varargout{2}=ve;
%varargout{3}=err;
%varargout{4}=count;


function[late,lone,R,V,u,v,ue,ve,err,errp,errn,count,flag,output]=eddyfit2d_one(r,cv,bool,z0,za,zb,lato,lono,options,str)

%vsize(r,cv,bool)

%count=sum(isfinite(cv),3);

cv=vmean(cv,3);
bool1=bool&isfinite(cv);
count=length(find(bool1));


if count>=4
    %vsize(cv,r)
    
    % %with angle from -pi to pi, angle guess = 0
    % [xf1,err1,flag1,output1]=fminsearchbnd( ...
    %     @(z) mean(squared(cv(bool1)-simpleddy(r(bool1)-z(2).*rot(z(1)),z(3),z(4),str))),...
    %     [0 z0],[-pi za],[pi zb],options);
    %
    % %with angle from 0 to 2pi, angle guess = pi
    % [xf2,err2,flag2,output2]=fminsearchbnd( ...
    %     @(z) mean(squared(cv(bool1)-simpleddy(r(bool1)-z(2).*rot(z(1)),z(3),z(4),str))),...
    %     [pi z0],[0 za],[2*pi zb],options);
    
    %with angle from -pi to pi, angle guess = 0
    [xf1,err1,flag1,output1]=fminsearchbnd( ...
        @(z) mean(squared(cv(bool1)-simpleddy(r(bool1)-z(2).*rot(z(1)),z(3),z(4),str))),...
        [0 z0],[-pi za],[pi zb],options);
    
    %with angle from 0 to 2pi, angle guess = pi
    [xf2,err2,flag2,output2]=fminsearchbnd( ...
        @(z) mean(squared(cv(bool1)-simpleddy(r(bool1)-z(2).*rot(z(1)),z(3),z(4),str))),...
        [pi z0],[0 za],[2*pi zb],options);
    
    if err1 < err2
        err=sqrt(err1);
        flag=flag1;
        output=output1;
        re=xf1(2).*rot(xf1(1));
        R=xf1(3);
        V=xf1(4);
    else
        err=sqrt(err2);
        flag=flag2;
        output=output2;
        re=xf2(2).*rot(xf2(1));
        R=xf2(3);
        V=xf2(4);
    end
else
   err=nan;
   flag=nan;
   output=nan;
   re=nan;
   R=nan;
   V=nan;
end

cve=simpleddy(r-re,R,V);
%cve(~bool)=nan;
ue=real(cve);
ve=imag(cve);

delta=(cv-cve).*rot(-angle(cve));
%vsize(cv,r,cve,bool1)

errp=sqrt(mean(squared(real(delta(bool1)))));
errn=sqrt(mean(squared(imag(delta(bool1)))));

%cv=vmean(cv,3);
u=real(cv);
v=imag(cv);
v(isnan(u))=nan;

[late,lone]=xy2latlon(real(re),imag(re),lato,lono);

%vsize(late,lone,R,V,ue,ve,err,count,flag,output)

function[]=eddyfit2d_test

z0=[0 80 -30];
za=[0 0 -80];
zb=[D inf 80];



%reporttest('EDDYFIT2D',aresame())

