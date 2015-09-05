function[bool]=inellipse(varargin)
%INELLIPSE  Locates points on the interior of ellipses.
%
%   INELLIPSE(Z,KAPPA,LAMBDA,THETA,ZO) where Z is a complex-valued
%   location Z=X+iY, and KAPPA, LAMBDA, THETA describe an ellipse at 
%   complex-valued position ZO=XO+iYO, is true if Z is inside the ellipse.
%
%   KAPPA is the ellipse amplitude, LAMBDA is the linearity, and THETA is
%   the orientation of the major axis with respect to the X-axis.
%
%   KAPPA and LAMBDA are related to the semi-axis lengths A and B by 
%   KAPPA^2=(A^2+B^2)/2 and LAMBDA=(A^2-B^2)/(A^2+B^2).  For details, see 
%   Lilly and Gascard (2006). 
%   __________________________________________________________________
%
%   Array or matrix input
%
%   BOOL=INELLIPSE(Z,KAPPA,LAMBDA,THETA,ZO), where Z is a length M *row*
%   vector, and the other input variables are length N *row* vector, 
%   returns an M x N boolean matrix BOOL.  
% 
%   In this case the (m,n)th entry of BOOL is true if the mth location Z
%   is inside the mth ellipse.
%
%   If the input variables have the same number of multiple *rows*, say K, 
%   this is interpreted as corresponding to K different time.  In this case
%   BOOL will be a 3D array with dimensions M x N x K.  
%   __________________________________________________________________
%
%   See also ELLCURVES, ELLIPSEPLOT.
%
%   'inellipse --f' generates sample figure.
%
%   Usage: bool=inellipse(z,kappa,lambda,theta,zo);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--f')
    type makefigs_inellipse
    makefigs_inellipse;
    return
end
 
z=varargin{1};
kappa=varargin{2};
lambda=varargin{3};
theta=varargin{4};
zo=varargin{5};

sizeall=vsize(z,zo,kappa,lambda,theta);
if ~allall((sizeall(:,1)-size(z,1))==0)
    error('All input arguments should have the same number of rows.')
end

bool=false(size(z,2),size(kappa,2),size(z,1));
for i=1:size(z,1)
    bool(:,:,i)=inellipse_one(z(i,:),zo(i,:),kappa(i,:),lambda(i,:),theta(i,:));
end

function[bool]=inellipse_one(z,zo,kappa,lambda,theta)

z=z(:);

%First we do an initial sort... not necessary
%dx=abs(a.*cos(theta).*cos(phi)-b.*sin(theta).*sin(phi));
%dy=abs(a.*sin(theta).*cos(phi)+b.*cos(theta).*sin(phi));

M=length(z);
N=length(kappa);

%Make all variables into matrices of the same size
z=vrep(z,N,2);
[zo,kappa,lambda,theta]=vrep(zo,kappa,lambda,theta,M,1);
[a,b]=kl2ab(kappa,lambda);


%Position relative to all ellipse centers
z=(z-zo);
vartheta=angle(z);


%Convert particle angle to ellipse phase
phi=atan(frac(a,b).*tan(vartheta-theta));

%Distance at this ellipse phase
R=kappa.*sqrt(1+lambda.*cos(2*phi));

bool=(abs(z)<=R);

%Find angle to ellipse center
