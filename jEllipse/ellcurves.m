function[x,y]=ellcurves(varargin)
%ELLCURVES  Returns curves corresponding to specified ellipse properties.
%
%   For direct plotting of ellipses, see ELLPLOT, which calls ELLCURVES.
%
%   ZC=ELLCURVES(KAPPA,LAMBDA,THETA,ZRES) returns complex-valued curves ZC
%   tracing out the periphery of ellipses with amplitude KAPPA, linearity 
%   LAMBDA, and orientation THETA, located at complex positions ZRES.
%
%   All input arguments are arrays of the same length, say N.  By default,  
%   ZC is calculated at 32 locations around the periphery.  Thus ZC will be
%   a complex-valued matrix with 32 rows and N columns. 
%
%   [XC,YC]=ELLCURVES(...) alternately returns the real and imaginary parts
%   of ZC, corresponding to the X- and Y- components of the ellipses.
%
%   ELLCURVES(KAPPA,LAMBDA,THETA,PHI,ZRES) with five input arguments begins 
%   each ellipse at phase PHI.  For most applications, the starting phase
%   does not matter.  The default value is to begin with a phase of zero. 
%
%   ELLCURVES(... ,'npoints',M) uses M points around the ellipse periphery.
%
%   ELLCURVES(... ,'aspect',AR) with AR=[XAR YAR] multiplies the X-signal 
%   and Y-signal by XAR and YAR, respectively, for plotting purposes. 
%
%   ELLCURVES is called by ELLPLOT, and is also useful in applications,
%   e.g. analysis of model fields within elliptical contours.
%
%   The curves can be plotted with PLOT(ZC).
%
%   See also ELLIPSEPLOT, ELLSIG, INELLIPSE.
%
%   Usage: zc=ellcurves(kappa,lambda,theta,zres);
%          zc=ellcurves(kappa,lambda,theta,phi,zres);
%          zc=ellcurves(kappa,lambda,theta,phi,zres,'npoints',64);
%          [xc,yc]=ellcurves(kappa,lambda,theta,phi,zres,'npoints',64);
%          zc=ellcurves(kappa,lambda,theta,phi,zres,'aspect',ar,'npoints',64);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2018 J.M. Lilly --- type 'help jlab_license' for details
 
 
na=length(varargin);

k=varargin{1};
l=varargin{2};
theta=varargin{3};

ar=[1 1];
npoints=32;

if na>4
    for i=1:3
        if ischar(varargin{end-1})
            if strcmpi(varargin{end-1}(1:3),'npo')
                npoints=varargin{end};
            elseif strcmpi(varargin{end-1}(1:3),'asp')
                ar=varargin{end};
            end
            na=na-2;
            varargin=varargin(1:end-2);
        end
    end
end

if length(varargin)==5
    phi=varargin{4};
    x=varargin{5};
else
    x=varargin{4};
    if iscell(k)
        for i=1:length(k)
             phi{i}=zeros(size(k{i}));
        end
    else
        phi=zeros(size(k));
    end 
end

if ~iscell(k)
    z=ellcurves_one(k,l,theta,phi,x,ar,npoints);
else
    for i=1:length(k)
        z{i,1}=ellcurves_one(k{i},l{i},theta{i},phi{i},x{i},ar,npoints);
    end
end


if nargout==2
    if iscell(z)
        x=cellreal(z);
        y=cellimag(z);
    else
        x=real(z);
        y=imag(z);
    end
else
    x=z;
end


function[z]=ellcurves_one(kappa,lambda,theta,phi,x,ar,npoints)
phio=linspace(0,2*pi,npoints);
%x=conj(x');

z=zeros(length(phio),length(kappa));
for i=1:length(phio)  %Loop over phio, very fast
    [xc,yc]=ellsig(kappa,lambda,theta,phi+phio(i),'real');
    %vsize(xc,yc,x)
    z(i,:)=ar(1)*xc+sqrt(-1)*ar(2)*yc+x;
end
%vsize(xc,yc,x)

%if ~isempty(strfind(str,'not'))
%    z=permute(z,[2 1]);  %This is like 6 times faster than conjugate transpose
%end
    

