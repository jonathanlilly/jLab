function[kappa,lambda,theta,phi,alpha,beta]=ellparams(varargin)
%ELLPARAMS  Ellipse parameters of a modulated bivariate or trivariate oscillation.
%
%   [KAPPA,LAMBDA,THETA,PHI]=ELLPARAMS(X,Y) where X and Y are analytic 
%   signals, returns the parameters of the complex-valued signal 
%   Z=REAL(X)+i REAL(Y), expressed as a modulated ellipse.
%
%   Here KAPPA is the RMS ellipse amplitude, LAMBDA is the linearity, 
%   THETA is the orientation, and PHI is the instantaneous orbital phase.
%
%   ELLPARAMS(M), where M is matrix with two columns, also works.
%
%   See Lilly and Gascard (2006) and Lilly and Olhede (2010a) for details.
%
%   ELLPARAMS is inverted by ELLSIG, which returns the X and Y signals 
%   given the ellipse parameters.
%
%   ELLPARAMS(...,DIM) performs the analysis with time running along
%   dimension DIM, as opposed to the default behavior of DIM=1.    
%   _______________________________________________________________________
%
%   Trivariate signals
%
%   ELLPARAMS also works for trivariate signals, which can be expressed as
%   a modulated ellipse in three dimensions.
%
%   [KAPPA,LAMBDA,THETA,PHI,ALPHA,BETA]=ELLPARAMS(X,Y,Z), where X, Y, and Z
%   are all analytic signals, also returns the zenith angle ALPHA and the 
%   azimuth angle BETA in addition to the other ellipse parameters.
%
%   ELLPARAMS(M), where M is matrix with three columns, also works.
%
%   See Lilly (2010) for details on the trivariate case.
%   _______________________________________________________________________
%
%   Cell array input / output
%
%   [KAPPA,LAMBDA,THETA,PHI]=ELLPARAMS(C) also works if C is a cell array
%   containing, say K different X and Y signals, each as a 2-column matrix,
%
%         C{1}(:,1)=X1, C{1}(:,2)=Y1
%         C{2}(:,1)=X2, C{2}(:,2)=Y2, ...  
%         C{K}(:,1)=XK, C{2}(:,2)=YK.  
%
%   In this case, the output variables will also be length K cell arrays.
%
%   The trivariate form described above also works, with each cell now
%   being a matrix with three columns.
%
%   This format works with the cell array output format of RIDGEWALK.
%   __________________________________________________________________
%
%   'ellparams --t' runs a test.
%
%   See also ELLSIG, ELLBAND, ELLDIFF, ELLVEL, ELLRAD, KL2AB, AB2KL. 
%
%   Usage: [kappa,lambda,theta,phi]=ellparams(x,y);
%          [kappa,lambda,theta,phi]=ellparams(x,y,dim);
%          [kappa,lambda,theta,phi,alpha,beta]=ellparams(x,y,z);
%          [kappa,lambda,theta,phi,alpha,beta]=ellparams([x,y,z]);
%          [kappa,lambda,theta,phi]=ellparams(C);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2016 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--t')
    ellparams_test,normvect_test,return
end

[na,k,l,theta,phi,alpha,beta]=vempty;

if ~iscell(varargin{end})&&length(varargin{end})==1
    dim=varargin{end};
    varargin=varargin(1:end-1);
else
    dim=1;
end

[z,kappa,lambda,theta,phi,alpha,beta]=vempty;
if ~isempty(varargin{1})
    if ~iscell(varargin{1})
        if length(varargin)==1
            if  (size(varargin{1},2)~=2)&&(size(varargin{1},2)~=3)
                error('For ELLPARAMS(M) with one input argument, M must have either two or three columns.')
            end
            x=varargin{1}(:,1);
            y=varargin{1}(:,2);
            if size(varargin{1},2)==3
                z=varargin{1}(:,3);
            end
        else
            if length(varargin)>3
                error('ELLPARAMS must have one, two, or three input arguments.')
            end
            x=varargin{1};
            y=varargin{2};
            if length(varargin)==3
                z=varargin{3};
            end
        end
        [kappa,lambda,theta,phi,alpha,beta]=ellparams_one(x,y,z,dim);
    else
        x=varargin{1};
        for i=1:length(x)
            if ~isempty(x{i})
                if size(x{i},2)==2
                    [kappa{i,1},lambda{i,1},theta{i,1},phi{i,1}]=ellparams_one(x{i}(:,1),x{i}(:,2),[],dim);
                else
                    [kappa{i,1},lambda{i,1},theta{i,1},phi{i,1},alpha{i,1},beta{i,1}]=...
                        ellparams_one(x{i}(:,1),x{i}(:,2),x{i}(:,3),dim);
                end
            else
                [kappa{i,1},lambda{i,1},theta{i,1},phi{i,1},alpha{i,1},beta{i,1}]=vempty;
            end
        end
    end
end

function[kappa,lambda,theta,phi,alpha,beta]=ellparams_one(x,y,z,dim)

[alpha,beta]=vempty;

if isempty(z)
    [kappa,lambda,theta,phi]=ellconv_xy2kl(abs(x),abs(y),angle(x),angle(y),dim);
else
    [nx,ny,nz]=normvect(x,y,z);
    warning('off','MATLAB:log:logOfZero')
    alpha=imag(log(sqrt(-1)*nx-ny));
    warning('on','MATLAB:log:logOfZero')
    beta=imag(log(nz+sqrt(-1)*sqrt(nx.^2+ny.^2)));
    [x,y,z]=vectmult(jmat3(-alpha,3),x,y,z);
    [x,y,z]=vectmult(jmat3(-beta,1),x,y,z);
    [kappa,lambda,theta,phi]=ellconv_xy2kl(abs(x),abs(y),angle(x),angle(y),dim);
end
      
function[kappa,lambda,theta,phi]=ellconv_xy2kl(X,Y,phix,phiy,dim)
%phia=double(phix+phiy+pi/2)/2;
%phid=double(phix-phiy-pi/2)/2;
%P=double(frac(1,2)*sqrt(squared(X)+squared(Y)+2.*X.*Y.*cos(2*phid)));
%N=double(frac(1,2)*sqrt(squared(X)+squared(Y)-2.*X.*Y.*cos(2*phid)));

phia=(phix+phiy+pi/2)/2;
phid=(phix-phiy-pi/2)/2;

P=frac(1,2)*sqrt(squared(X)+squared(Y)+2.*X.*Y.*cos(2*phid));
N=frac(1,2)*sqrt(squared(X)+squared(Y)-2.*X.*Y.*cos(2*phid));

phip=unwrap(phia+imlog(X.*rot(phid)+Y.*rot(-phid)),[],dim);
phin=unwrap(phia+imlog(X.*rot(phid)-Y.*rot(-phid)),[],dim);

kappa=sqrt(P.^2+N.^2);
lambda=frac(2*P.*N.*sign(P-N),P.^2+N.^2);

%For vanishing linearity, put in very small number to have sign information 
lambda(lambda==0)=sign(P(lambda==0)-N(lambda==0))*(1e-10);

theta=phip/2-phin/2;
phi=  phip/2+phin/2;

theta=unwrap(theta,[],dim);
phi=unwrap(phi,[],dim);

lambda=real(lambda);


function[nx,ny,nz]=normvect(x,y,z)
%NORMVECT  Unit normal vector to the ellipse plane in three dimensions.
%
%   [NX,NY,NZ]=NORMVECT(X,Y,Z) returns the three components of the unit 
%   normal vector to the plane containing the ellipse specified by the 
%   three complex-valued arrays X, Y, and Z.
%
%   In vector notation the normal vector is defined as N=IMAG(X) x REAL(X),
%   where ``x'' is the vector cross product, and the unit normal is N/||N||.
%
%   The input arrays and output arrays are all the same size.
% 
%   See Lilly (2010) for details.
%
%   'normvect --t' runs a test.
%
%   Usage: [nx,ny,nz]=normvect(x,y,z);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2010 J.M. Lilly --- type 'help jlab_license' for details
 
nx= (imag(y).*real(z)-imag(z).*real(y));
ny=-(imag(x).*real(z)-imag(z).*real(x));
nz= (imag(x).*real(y)-imag(y).*real(x));
  
denom=sqrt(nx.^2+ny.^2+nz.^2);
nx=frac(nx,denom);
ny=frac(ny,denom);
nz=frac(nz,denom);

function[]=normvect_test
 
load solomon 
use solomon

[x,y,z]=anatrans(x./1e4,y./1e4,z./1e4);
[nx,ny,nz]=normvect(x,y,z);

dot=nx.*x+ny.*y+nz.*z;

reporttest('NORMVECT parallel part of X_+ vanishes, Solomon Islands',allall(abs(dot)<1e-10))


%Choose central part where signal is elliptical
vindex(x,y,z,7000:10000,1);

[ax,omx,upx]=instmom(x);
[ay,omy,upy]=instmom(y);
[az,omz,upz]=instmom(z);
     
dx=x.*(upx+sqrt(-1)*omx);
dy=y.*(upy+sqrt(-1)*omy);
dz=z.*(upz+sqrt(-1)*omz);

[nx,ny,nz]=normvect(x,y,z);

dot=nx.*dx+ny.*dy+nz.*dz; %Parallel part of derivative

[dnx,dny,dnz]=vdiff(nx,ny,nz,1);

dot2=-(dnx.*x+dny.*y+dnz.*z); 

err=abs(dot-dot2).^2./abs(dot).^2;
err=flipud(sort(err));
err=err(60:end);


reporttest('NORMVECT derivative of parallel part matches, Solomon Islands (removing worst outliers)',allall(err<0.05))




function[]=ellparams_test
 
t=(0:1:925)';
kappa=3*exp(2*0.393*(t/1000-1));
lambda=0.5+0*kappa;
phi=(t/1000*5)*2*pi;
theta=pi/4+phi./14.45;

beta=pi/6+phi./14.45*lambda(1);
alpha=pi/6-phi./14.45*2*lambda(1)*sqrt(2);

[x,y,z]=ellsig(kappa,lambda,theta,phi,alpha,beta);
[kappa2,lambda2,theta2,phi2,alpha2,beta2]=ellparams(x,y,z);

x1=[kappa lambda theta phi alpha beta];
x2=[kappa2 lambda2 theta2 phi2 alpha2 beta2];
reporttest('ELLPARAMS rapidly changing trivariate ellipse',aresame(x1,x2,1e-8))



% %/*******************************************************************
% %Flip ellipse parameters if unit normal is pointing radially inwards
% %Components of unit normal vector to surface of earth
% [nx,ny,nz]=normvect(xr,yr,zr);
% 
% %Replicate x, y, and z along columns
% xmat=vrep(x,size(nx,2),2)./radearth;
% ymat=vrep(y,size(ny,2),2)./radearth;
% zmat=vrep(z,size(nz,2),2)./radearth;
% 
% %Projection of normal vector to plane onto normal to sphere
% proj=xmat.*nx+ymat.*ny+zmat.*nz;
% 
% bool=(proj<0);
% theta(bool)=-theta(bool);
% lambda(bool)=-lambda(bool);
% beta(bool)=pi+beta(bool);
% nx(bool)=-nx(bool);
% ny(bool)=-ny(bool);
% nz(bool)=-nz(bool);
% 
% if length(find(bool))>0
%     [xr2,yr2,zr2]=ellsig(kappa,lambda,theta,phi,alpha,beta);
%     tol=1e-6;
%     reporttest('ELLIPSEXTRACT adjustment for sign of normal vector',aresame(xr2,xr,tol)&&aresame(yr2,yr,tol)&&aresame(zr2,zr,tol))
% end
% 
% %dev=1-abs(proj);
% %figure,plot(dev)
% [latn,lonn]=xyz2latlon(nx*radearth,ny*radearth,nz*radearth);
% %\*******************************************************************
