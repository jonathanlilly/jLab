function[varargout]=polysmooth(varargin)
%POLYSMOOTH  Smoothing scattered 2D data with local polynomial fitting.
%
%   POLYSMOOTH generates a map from scattered data in two dimensions---
%   either on the plane or on the sphere---using a local least squares fit 
%   to a polynomial. 
%
%   A numerically efficient algorithm is used that avoids explicit loops. 
%   Also, the data are pre-sorted so that different mapping parameters can 
%   be tried out at little computational expense.
%
%   A reference paper for this function is currently in preparation.  In 
%   the meantime, for more details see the following book:
%
%      Fan and Gijbels, 1996.  Local polynomial modelling and its 
%          applications. Chapman and Hall.
%   __________________________________________________________________
%
%   Smoothing on the plane
%
%   Let's say we have an array Z of data is at locations X,Y, where X,Y, 
%   and Z are all arrays of the same size.  The problem is to obtain a 
%   mapped field ZHAT on some regular grid specified by arrays XO and YO.
%
%   Calling POLYSMOOTH is a two-step process: 
%
%       [DS,XS,YS,ZS]=TWODSORT(X,Y,Z,XO,YO,CUTOFF);
%       ZHAT=POLYSMOOTH(DS,XS,YS,ZS,B,P);
%
%   Firstly, TWODSORT which returns ZS, a 3D array of data values at each 
%   grid point, sorted by increasing distance DS, and the corresponding 
%   positions XS and YS. Here XO and YO are arrays specifying the bin 
%   center locations of the grid.  See TWODSORT for more details.
%  
%   Here CUTOFF determines the maximum distance included in the sorting 
%   and should be chosen to be greater than B.  
%
%   Secondly, POLYSMOOTH fits a Pth order polynomial at each gridpoint 
%   within a neighborhood specified by the "bandwidth" B.  
%
%   Data and grid point locations are presumed to have the same units as
%   the bandwidth B (e.g., kilometers).
%
%   P may be chosen as P=0 (fit to a constant), P=1 (fit to a plane), or
%   else P=2 (fit to a parabolic surface).  
%
%   POLYSMOOTH can also be used for data on the sphere, as described next.
%   __________________________________________________________________
%
%   Smoothing on the sphere
%
%   POLYSMOOTH supports a local polynomial fit on the sphere.  As before 
%   this is a two-step process:
%
%       [DS,XS,YS,ZS]=SPHERESORT(LAT,LON,Z,LATO,LONO,CUTOFF);
%       ZHAT=POLYSMOOTH(DS,XS,YS,ZS,B,P,'sphere') 
%
%   Firstly one calls SPHERESORT, which is the analogue of TWODSORT for the 
%   sphere.  LATO and LONO are arrays specifying the bin center locations 
%   of the grid.  See SPHERESORT for more details.
%
%   Secondly, POLYSMOOTH fits a Pth order polynomial at each gridpoint 
%   within a neighborhood specified by the "bandwidth" B.  The bandwidth 
%   here should have units of kilometers. 
%
%   Note that SPHERESORT and POLYSMOOTH both assume the sphere to be the 
%   radius of the earth, as specified by the function RADEARTH.
%   __________________________________________________________________
%  
%   Additional options
%   
%   A number of different options are possible, as descibed below.  These  
%   are specified with trailing string arguments, which can be given in any
%   order provided they are after the numeric arguments.
%   __________________________________________________________________
%
%   Weighting function
%
%   POLYSMOOTH weights the data points in the vicinity of each grid point
%   according to some decaying function of distance.
%
%   1) Parabolic weighting
%
%   POLYSMOOTH(...,'epa'), the default behavior, uses the Epanechnikov  
%   kernel --- a fancy name for a parabolic weighting function.  This  
%   function vanishes outside a radius of the bandwidth B.  
%
%   The Epanechnikov kernel has the numerical advantage that distances from 
%   far away points do not have to be computed, and is therefore preferred.
%
%   2) Gaussian weighting:
%
%   POLYSMOOTH(...,'gau') uses a Gaussian weighting function instead, with 
%   a standard deviation of B and truncated, for numerical reasons, at a 
%   radius of 3B.  
%
%   When using the Gaussian weighting function you'll need to call TWODSORT
%   or SPHERESORT with a cutoff distance of at least 3B.
%   __________________________________________________________________
%
%   Algorithms 
%
%   POLYSMOOTH has two algorithm choices for optimal performance, both of 
%   which give identical answers.
%
%   1) Speed optimization
%
%   The default algorithm is optimized for speed but uses a great deal of 
%   memory.  It solves the least squares problem for all grid points
%   simultaneously (by directly solving matrix inversions---see MATINV.)
%
%   This can be greater than an order of magnitude faster than the obvious 
%   approach of looping over the grid points.
% 
%   2) Memory optimization
%
%   If the data is so large that memory becomes a limiting factor, use
%
%       ZHAT=POLYSMOOTH(DS,XS,YS,ZS,B,P,'memory');
%  
%   which performs an explicit loop. 
%   __________________________________________________________________
%
%   Weighted data points
%
%   POLYSMOOTH can also incorporate a weighting factor on the data points.
%   Let W be an array of positive values the same size as the data array Z.  
%   Each data point is then treated as if it were W data points.
%
%   To use weighted data points, call TWODSORT or SPHERESORT with W added:
%
%       [DS,XS,YS,ZS,WS]=TWODSORT(X,Y,Z,W,XO,YO,CUTOFF);           
%    or [DS,XS,YS,ZS,WS]=SPHERESORT(LAT,LON,Z,W,LATO,LONO,CUTOFF); 
%    
%   followed by    ZHAT=POLYSMOOTH(DS,XS,YS,ZS,WS,B,P);
%            or    ZHAT=POLYSMOOTH(DS,XS,YS,ZS,WS,B,P,'sphere');
%
%   for data distributed on the plane or on the sphere, respectively.
%
%   For very large datasets, this is useful in first condensing the data
%   into a manageable size.
%   __________________________________________________________________
%
%   Additional output arguments
%
%   [ZHAT,WEIGHT]=POLYSMOOTH(...) also returns the total weight used in 
%   the computation at each grid point.  WEIGHT is the same size as ZHAT.
%
%   Small weights mean that only few or far away data points contributed
%   to the calculation; therefore, the values of ZHAT at these points
%   are less reliable than at points where the total weight is large.
%
%   [ZHAT,WEIGHT,BETA]=POLYSMOOTH(...) additionally returns the matrix
%   of the estimated first P derivatives.  BETA is 3D and is the same 
%   size as ZHAT and WEIGHT along its first two dimensions.
%
%   BETA(:,:,1) is the estimated zeroth derivative, i.e. the surface, so
%   this is the same as ZHAT.
%   __________________________________________________________________
%
%   Variable bandwidth
%
%   ZHAT=POLYSMOOTH(...,N,P,'variable') uses the parameter N instead of a
%   fixed bandwidth.  The bandwidth now varies over the grid such that N
%   points fall within one bandwidth distance from every grid point.
%   
%   [ZHAT,WEIGHT,BETA,B]=POLYSMOOTH(...,N,P,'variable') returns the 
%   bandwidth B at each grid point, which is the same size as ZHAT.  
%
%   The variable bandwidth algorithm can give good results when the data
%   spacing is highly uneven, particularly when used with a higher-order
%   (P=1 or P=2) fit.  
%
%   When using the variable bandwidth method, be aware that the bandwidth
%   necessary to include N points may turn out to be larger than CUTOFF,
%   the maximum distance value given to TWODSORT or SPHERESORT.  
%   __________________________________________________________________
%
%   'polysmooth --t' runs some tests.
%   'polysmooth --f' generates some sample figures.
%   'polysmooth --f2' with jData installed generates the figure shown
%       above, which may require a relatively powerful computer.
%
%   Usage:  [ds,xs,ys,zs]=twodsort(x,y,z,xo,yo,cutoff);  
%           [zhat,weight,beta]=polysmooth(ds,xs,ys,zs,b,p);
%   --or--
%           [ds,xs,ys,zs]=spheresort(lat,lon,z,w,lato,lono,cutoff); 
%           [zhat,weight,beta]=polysmooth(ds,xs,ys,zs,b,p,'sphere');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    polysmooth_test,return
elseif strcmpi(varargin{1}, '--f')
    type makefigs_polysmooth
    makefigs_polysmooth;
    return
elseif strcmpi(varargin{1}, '--f2')
    type makefigs_polysmooth2
    makefigs_polysmooth2;
    return
end

str='speed';
kernstr='epanechnicov';
geostr='cartesian';
varstr='constant';

for i=1:4
    if ischar(varargin{end})
        tempstr=varargin{end};
        if strcmpi(tempstr(1:3),'epa')||strcmpi(tempstr(1:3),'gau')
            kernstr=tempstr;
        elseif strcmpi(tempstr(1:3),'car')||strcmpi(tempstr(1:3),'sph')
            geostr=tempstr;
        elseif strcmpi(tempstr(1:3),'spe')||strcmpi(tempstr(1:3),'loo')||strcmpi(tempstr(1:3),'mem')
            str=tempstr;
        elseif strcmpi(tempstr(1:3),'var')||strcmpi(tempstr(1:3),'con')
            varstr=tempstr;
        end
        varargin=varargin(1:end-1);
    end
end

R=varargin{end};
B=varargin{end-1};

na=length(varargin)-2;

if strcmpi(str(1:3),'spe')||strcmpi(str(1:3),'mem')
    d=varargin{1};
    x=varargin{2};
    y=varargin{3};
    z=varargin{4};
    if na==5
        w=varargin{5};
    else
        w=ones(size(z));
    end
elseif strcmpi(str(1:3),'loo')
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    if na==5
        w=ones(size(z));
        xo=varargin{4};
        yo=varargin{5};
    elseif na==6
        w=varargin{4};
        xo=varargin{5};
        yo=varargin{6};
    end
end

if strcmpi(str(1:3),'loo')
    [beta,weight,B]=polysmooth_slow(x,y,z,w,xo,yo,B,R,varstr,kernstr,geostr);
elseif strcmpi(str(1:3),'spe')
    [beta,weight,B]=polysmooth_fast(d,x,y,z,w,B,R,varstr,kernstr);
elseif strcmpi(str(1:3),'mem')
    [beta,weight,B]=polysmooth_memory(d,x,y,z,w,B,R,varstr,kernstr);
else
    error(['Algorithm type ' str ' is not supported.'])
end

%Adjustment for complex-valued data
if ~allall(isreal(z))
    beta(isnan(real(beta(:))))=nan+sqrt(-1)*nan;
end

varargout{1}=beta(:,:,1);
varargout{2}=weight;
varargout{3}=beta;
varargout{4}=B;


%for k=1:max(size(beta,3),nargout)
%    varargout{k}=squeeze(beta(:,:,k));
%end

function[B]=polysmooth_bandwidth(varstr,B,M,N,d,w)

if strcmpi(varstr(1:3),'con')
    if size(B,2)==1
        B=vrep(B,N,2);
    end
    if size(B,1)==1
        B=vrep(B,M,1);
    end
else
%%Earlier version before weighting
%     if ndims(d)<3
%         B=d(B);
%     else
%         B=squeeze(d(:,:,B));
%     end

     N=B;
     %Compute number of data points with potential for weighted points
     if  ndims(d)<3
         B=nan;
         cumw=cumsum(double(w),1);
         index=find(cumw>=N,1,'first');
         if ~isempty(index)
              B=squeeze(d(index));
         end
     else
         B=nan*zeros(size(d,1),size(d,2));
         cumw=cumsum(double(w),3);
         for i=1:size(d,1);
             for j=1:size(d,2)
                 index=find(cumw(i,j,:)>=N,1,'first');
                 if ~isempty(index)
                      B(i,j)=squeeze(d(i,j,index));
                 end
             end
         end
     end
end

function[beta,weight,B]=polysmooth_fast(ds,xs,ys,zs,ws,B,R,varstr,kernstr)
tol=1e15;

M=size(xs,1);
N=size(xs,2);

B=polysmooth_bandwidth(varstr,B,M,N,ds,ws);
B=vrep(B,size(xs,3),3);

P=sum(0:R+1);

matcond=zeros(M,N);


beta=nan*zeros(size(xs,1),size(xs,2),P);
weight=nan*zeros(size(xs,1),size(xs,2));
%Initial search to exclude right away points too far away
ds(isnan(zs)|ds>B)=nan;
numgood=squeeze(sum(sum(~isnan(ds),2),1));
index=find(numgood~=0,1,'last'); 

if ~isempty(index)
    vindex(ds,xs,ys,zs,ws,B,1:index,3);    
    W=polysmooth_kern_dist(ds,B,kernstr).*ws;  %ws for weighted data ponts
    vswap(ds,xs,ys,zs,W,nan,0);
    weight=sum(W,3);
    
    X=polysmooth_xmat(xs,ys,R);
 
    XtimesW=X.*vrep(W,P,4);
    mat=zeros(M,N,P,P);
    vect=zeros(M,N,P);
    for i=1:size(X,4)
        mat(:,:,:,i)=squeeze(sum(XtimesW.*vrep(X(:,:,:,i),P,4),3));
    end
    vect=squeeze(sum(XtimesW.*vrep(zs,P,4),3));
  
    for i=1:N
         for j=1:M
                 matcond(j,i)=cond(squeeze(mat(j,i,:,:)));
                 if R==1
                    if matcond(j,i).^2>tol
                           mat(j,i,:,:)=nan;
                    end
                 elseif R==2
                    if matcond(j,i)>tol
                           mat(j,i,:,:)=nan;
                    end
                 end
         end
    end
 
    if P==1
        vswap(mat,0,nan);
        beta=vect./mat;
    else
        
        beta=matmult(matinv(mat),vect,3);
                
%         sizemat=size(mat);
%         mat=reshape(mat,size(mat,1).*size(mat,2),size(mat,3),size(mat,4));
%         vect=reshape(vect,size(vect,1).*size(vect,2),size(vect,3));
%         nonnani=~isnan(sum(sum(mat,2),3));
%         beta=nan*zeros(size(mat,1),size(mat,2)); 
%         
%         vindex(mat,vect,nonnani,1); 
%         size(mat)
%         beta(nonnani,:)=matmult(matinv(mat),vect,2);
%         beta=reshape(beta,sizemat(1),sizemat(2),sizemat(3));
%         vsize(beta,beta1)
%         
%         aresame(beta1,beta,1e-6)
%         toc
    end    
end
weight=matcond;
B=B(:,:,1);

function[beta,weight,B]=polysmooth_memory(d,x,y,z,w,B,R,varstr,kernstr)

tol=1e-10;
P=sum(0:R+1);

M=size(x,1);
N=size(x,2);
    
B=polysmooth_bandwidth(varstr,B,M,N,d,w);

beta=nan*zeros(M,N,P);
weight=nan*zeros(M,N,2);
for j=1:M
    disp(['Polysmooth computing map for row ' int2str(j) '.'])
    for i=1:N
        [dist,xd,yd,zd,wd]=vsqueeze(d(j,i,:),x(j,i,:),y(j,i,:),z(j,i,:),w(j,i,:));
        
        W=polysmooth_kern_dist(dist,B(j,i),kernstr);
       
        index=find(W~=0&~isnan(zd));       
        if ~isempty(index)
            vindex(xd,yd,zd,wd,index,1);
            
            weight(j,i)=sum(W(index));
            W=diag(W(index).*wd);  %Times wd for weighted dat points
                
            X=squeeze(polysmooth_xmat(xd,yd,R));
            
            XW=conj(X')*W;
            mat=XW*X;
            matcond=rcond(mat);
            if matcond>tol;
                %beta(j,i,:)=inv(mat)*XW*zd;
                beta(j,i,:)=mat\(XW*zd);  %Following Checkcode's suggestion
            end
        end
    end
end

function[beta,weight,B]=polysmooth_slow(x,y,z,w,xo,yo,B0,R,varstr,kernstr,geostr)

%Explicit loop without pre-sorting, for testing purposes only
M=length(yo);
N=length(xo);
    
vsize(x,y,z,xo,yo,B0,R);

tol=1e-10;
P=sum(0:R+1);

beta=nan*zeros(M,N,P);
weight=nan*zeros(M,N);
B=nan*zeros(M,N);

for j=1:M
    %Note that the slow method loops over the other direction for sphere
    %than the memory method
    disp(['Polysmooth computing map for row ' int2str(j) '.'])
    for i=1:N
        if strcmpi(geostr(1:3),'sph')
            [xd,yd,dist]=latlon2xy(x,y,xo(i),yo(j));            
        else
            xd=x-xo(i);
            yd=y-yo(j);
            dist=sqrt(xd.^2+yd.^2);
        end
        zd=z;
        wd=w;
        
        [sortdist,sorter]=sort(dist);
        B(j,i)=polysmooth_bandwidth(varstr,B0,1,1,sortdist,wd(sorter));
        W=polysmooth_kern_dist(dist,B(j,i),kernstr);

        %size(W),size(zd)
        index=find(W~=0&~isnan(zd));       
        if ~isempty(index)
                  
            vindex(xd,yd,zd,wd,index,1);
            X=squeeze(polysmooth_xmat(xd,yd,R));
            
            weight(j,i)=sum(W(index));
            W=diag(W(index).*wd);  %Times wd for weighted data points
            
            XW=conj(X')*W;
            mat=XW*X;
            matcond=rcond(mat);
            if matcond>tol;
                %beta(j,i,:)=inv(mat)*XW*zd;
                beta(j,i,:)=mat\(XW*zd); %Following Checkcode's suggestion
            end
        end
    end
end

if strcmpi(geostr(1:3),'sph')
    beta=permute(beta,[2 1 3]);
    weight=permute(weight,[2 1]); 
    B=permute(B,[2 1]); 
end

function[X]=polysmooth_xmat(x,y,R)
P=sum(0:R+1);
X=ones(size(x,1),size(x,2),size(x,3),P);
if R>0
    X(:,:,:,2)=x;
    X(:,:,:,3)=y;
end
if R>1
    X(:,:,:,4)=frac(1,2)*x.^2;   %See writeup
    X(:,:,:,5)=x.*y;   
    X(:,:,:,6)=frac(1,2)*y.^2;   %See writeup
end


% function[W]=polysmooth_kern_xy(x,y,Binv,kernstr);
% 
% 
% if size(Binv,4)==1
%     xp=Binv.*x;
%     yp=Binv.*y;
% elseif size(Binv,4)==2
%     xp=Binv(:,:,:,1).*x;
%     yp=Binv(:,:,:,2).*y;
% %elseif size(Binv,3))==4
% %    xp=Binv(1,1).*x+Binv(1,2).*y;
%  %   yp=Binv(2,1).*x+Binv(2,2).*y;
% end
% 
% distsquared=abs(xp).^2+abs(yp).^2;
% 
% if strcmpi(kernstr(1:3),'epa')
%    W=frac(2,pi).*(1-distsquared).*(1+sign(1-sqrt(distsquared)));
% elseif strcmpi(lower(kernstr(1:3)),'gau')
%    W=frac(2,pi).*exp(-frac(distsquared,2));%Not right, should have frac(1,sigx*sigy)
% end

function[W]=polysmooth_kern_dist(dist,B,kernstr)

% if iscell(B)
%     N=B{2};
%     B=B{1};
% else
%     N=3;
% end
N=3;

dist=frac(dist,B);

if strcmpi(kernstr(1:3),'epa')
   W=frac(2,pi).*(1-dist.^2).*(1+sign(1-dist));
elseif strcmpi(kernstr(1:3),'gau')
   W=frac(2,pi.*B.^2).*exp(-frac(dist.^2,2)).*(1+sign(1-dist./N));
end

W=vswap(W,nan,0);
        
function[]=polysmooth_test
 
polysmooth_test_cartesian;
polysmooth_test_sphere;


function[]=polysmooth_test_cartesian
 
%Use peaks for testing... random assortment
[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:200);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.5:3);
yo=(-3:.6:3);

B=2;

[ds,xs,ys,zs]=twodsort(xdata,ydata,zdata,xo,yo,B);    
tic;z1=polysmooth(xdata,ydata,zdata,xo,yo,B,0,'epan','loop');etime1=toc;
tic;z2=polysmooth(ds,xs,ys,zs,B,0,'epan','memory');etime2=toc;
tic;z3=polysmooth(ds,xs,ys,zs,B,0,'epan','speed');etime3=toc;

disp(['POLYSMOOTH was ' num2str(etime1./etime3) ' times faster than loop, constant fit.'])
disp(['POLYSMOOTH was ' num2str(etime2./etime3) ' times faster than memory method, constant fit.'])

tol=1e-8;
b1=aresame(z1,z3,tol);
b2=aresame(z2,z3,tol);

reporttest('POLYSMOOTH output size matches input size',aresame(size(z1),[size(xs,1) size(xs,2)]))
reporttest('POLYSMOOTH speed and loop methods are identical for constant fit',b1)
reporttest('POLYSMOOTH speed and memory methods are identical for constant fit',b2)

tic;[z1,w1,beta1]=polysmooth(xdata,ydata,zdata,xo,yo,B,1,'epan','loop');etime1=toc;
tic;[z2,w2,beta2]=polysmooth(ds,xs,ys,zs,B,1,'epan','memory');etime2=toc;
tic;[z3,w3,beta3]=polysmooth(ds,xs,ys,zs,B,1,'epan','speed');etime3=toc;

disp(['POLYSMOOTH was ' num2str(etime1./etime3) ' times faster than loop, linear fit.'])
disp(['POLYSMOOTH was ' num2str(etime2./etime3) ' times faster than memory method, linear fit.'])

tol=1e-8;
b1=aresame(beta1,beta3,tol);
b2=aresame(beta2,beta3,tol);

reporttest('POLYSMOOTH speed and loop methods are identical for linear fit',b1)
reporttest('POLYSMOOTH speed and memory methods are identical for linear fit',b2)


function[]=polysmooth_test_sphere
 
%Use peaks for testing... random assortment
[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:200);

%Convert to Lat and Longitude
lon=x.*60;
lat=y.*30;
[latdata,londata,zdata]=vindex(lat(:),lon(:),z(:),index,1);
wdata=10*abs(rand(size(zdata)));  %For heavy data point test

lono=(-180:20:180);
lato=(-90:20:90);

B=2000;

[ds,xs,ys,zs,ws]=spheresort(latdata,londata,zdata,wdata,lato,lono,B,'loop');

tic;z1=polysmooth(latdata,londata,zdata,lato,lono,B,0,'epan','sphere','loop');etime1=toc;
tic;z2=polysmooth(ds,xs,ys,zs,B,0,'epan','memory');etime2=toc;
tic;z3=polysmooth(ds,xs,ys,zs,B,0,'epan','speed');etime3=toc;

disp(['POLYSMOOTH was ' num2str(etime1./etime3) ' times faster than loop, constant fit on a sphere.'])
disp(['POLYSMOOTH was ' num2str(etime2./etime3) ' times faster than memory method, constant fit on a sphere.'])

tol=1e-8;
b1=aresame(z1,z3,tol);
b2=aresame(z2,z3,tol);

reporttest('POLYSMOOTH output size matches input size on a sphere',aresame(size(z1),[size(xs,1) size(xs,2)]))
reporttest('POLYSMOOTH speed and loop methods are identical for constant fit on a sphere',b1)
reporttest('POLYSMOOTH speed and memory methods are identical for constant fit on a sphere',b2)

tic;z1=polysmooth(latdata,londata,zdata,wdata,lato,lono,B,0,'epan','sphere','loop');etime1=toc;
tic;z2=polysmooth(ds,xs,ys,zs,ws,B,0,'epan','memory');etime2=toc;
tic;z3=polysmooth(ds,xs,ys,zs,ws,B,0,'epan','speed');etime3=toc;

tol=1e-8;
b1=aresame(z1,z3,tol);
b2=aresame(z2,z3,tol);

reporttest('POLYSMOOTH speed and loop methods are identical with heavy data points',b1)
reporttest('POLYSMOOTH speed and memory methods are identical with heavy data points',b2)

B0=10;

[ds,xs,ys,zs]=spheresort(latdata,londata,zdata,lato,lono,radearth*pi/4-1e-3);

tic;[z1,weight1,beta1,B1]=polysmooth(latdata,londata,zdata,lato,lono,B0,0,'epan','sphere','variable','loop');etime1=toc;
tic;[z2,weight2,beta2,B2]=polysmooth(ds,xs,ys,zs,B0,0,'epan','memory','variable');etime2=toc;
tic;[z3,weight3,beta3,B3]=polysmooth(ds,xs,ys,zs,B0,0,'epan','speed','variable');etime3=toc;

disp(['POLYSMOOTH was ' num2str(etime1./etime3) ' times faster than loop, constant 10-point fit on a sphere.'])
disp(['POLYSMOOTH was ' num2str(etime2./etime3) ' times faster than memory method, constant 10-point fit on a sphere.'])

tol=1e-8;
%There are some NANs in z3 that do not appear in z1, for some reason
b1=aresame(z1(isfinite(z3)),z3(isfinite(z3)),tol);
b2=aresame(z2,z3,tol);
%figure,plot(z1,'b'),hold on,plot(z3,'r.')

reporttest('POLYSMOOTH output size matches input size on a sphere',aresame(size(z1),[size(xs,1) size(xs,2)]))
reporttest('POLYSMOOTH speed and loop methods are identical for constant 10-point fit on a sphere',b1)
reporttest('POLYSMOOTH speed and memory methods are identical for constant 10-point fit on a sphere',b2)



% %Zeroth-order fit at 200~km radius
% [zhat0,weight,beta,b]=polysmooth(ds,xs,ys,zs,200,0,'sphere');
% 
% figure
% contourf(lono,lato,zhat0,[0:1/2:50]),nocontours,caxis([0 45]),
% latratio(30),[h,h2]=secondaxes;topoplot continents,
% title('Standard Deviation of SSH from ALONGTRACK.MAT, mapped using 200~km smoothing')
% hc=colorbar('EastOutside');axes(hc);ylabel('SSH Standard Deviation (cm)')
% set(h,'position',get(h2,'position'))
% 
% currentdir=pwd;
% cd([whichdir('jlab_license') '/figures'])
% print -dpng alongtrack_std_constant
% crop alongtrack_std.png
% cd(currentdir)