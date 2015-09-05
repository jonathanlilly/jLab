function[varargout]=anatrans(varargin)
%ANATRANS  Analytic part of signal.
%
%   XP=ANATRANS(X) returns the analytic part of the real-valued signal X,
%   which is a column vector or a matrix with 'time' in columns. 
%
%   XP is defined for real X such that X = REAL(XP) = 1/2*(XP+CONJ(XP)).
% 
%   ANATRANS removes the mean of the signal.  It does not detrend, so if 
%   your time series has a linear trend you may wish to call DETREND first. 
%
%   [XP1,XP2,...,XPM]=ANATRANS(X1,X2,...,XM) also works for multiple input
%   arguments.  
%    
%   ANATRANS(X,DIM) optionally applies the analytic transform along 
%   dimension DIM, instead of along the first dimension.
%   ___________________________________________________________________
%   
%   Complex-valued signals
%
%   [ZP,ZN]=ANATRANS(Z,CONJ(Z)) returns the analytic part ZP of the 
%   complex-valued signal Z, and the analytic part ZN of its conjugate.
%
%   The normalization for complex-valued signals differs from that for 
%   real valued signals.  ZP and ZN are defined such that Z=ZP+CONJ(ZN). 
%
%   The analytic part of a real-valued signal and that of a complex-valued 
%   signal are thus defined differently by a factor of two, following the 
%   convention of Lilly and Gascard (2006) and Lilly and Olhede (2010).
%
%   Note that the equality Z=ZP+CONJ(ZN) is not necessarily exact, due to 
%   issues relating to boundary effects from the time series edges.
%   ___________________________________________________________________
%
%   Boundary conditions
%
%   ANATRANS(...,STR), where STR is a string, optionally specifies the
%   boundary condition to be imposed at the edges of the time series.  
%   Valid options for STR are 
%
%         STR = 'periodic' for periodic boundary conditions 
%         STR = 'zeros' for zero-padding beyond the endpoints 
%         STR = 'mirror' for reflecting the time series at both ends
%
%   The default value of STR is 'mirror', as this tends to minimize the
%   `edge effects' that occur near the ends of the time series.
%   ___________________________________________________________________
%
%   'anatrans --t' runs some tests.
%
%   Usage: z=anatrans(x);
%          z=anatrans(x,dim);
%          z=anatrans(x,dim,'mirror');
%          [zp,zn]=anatrans(z,conj(z));
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        


if strcmpi(varargin{1},'--t')
    anatrans_test,return
end

% antrans is basically the same as Matlab's Hilbert:
% if isreal(x)
%     zp=x+sqrt(-1)*hiltrans(x);
% else
%     zp=frac(1,2)*(x+sqrt(-1)*hiltrans(x));
%     zn=frac(1,2)*(x-sqrt(-1)*hiltrans(x));
% end

str='mirror';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end
dim=1;
if length(varargin{end})==1
    dim=varargin{end};
    varargin=varargin(1:end-1);
end

for i=1:length(varargin)
    varargout{i}=anatrans_one(varargin{i},str,dim);
end

function[zp]=anatrans_one(x,str,dim)

M0=size(x,dim);

mx=vmean(x,dim);
x=x-vrep(mx,M0,dim);
x=timeseries_boundary(x,dim,str,'nodetrend');
M=size(x,dim);

Z=2*fft(x,[],dim);

index=round(M/2):size(Z,dim);

Z=vindexinto(Z,0,index,dim);
zp=ifft(Z,[],dim);

if ~isreal(x)
    zp=frac(1,2)*zp;
end

if M0~=M
    index=M0+1:M0*2;
    vindex(zp,index,dim);     
end

function[y]=hiltrans(x)
%HILTRANS  Hilbert transform.
%
%   Y=HILTRANS(X) returns the Hilbert transform of column vector X.
%
%   If SIZE(X,2)>1, HILTRANS takes the Hilbert transform along columns.
%
%   'hiltrans --f' makes a sample figure
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        
  
if isreal(x)
  bool=1;
else
  bool=0;
end

N=size(x,1);
X=fft(x);
sgnom=ones(N,1);
%sgnom(1)=0;
index=(1:N);
sgnom(index-1>N/2)=-1;
d=-sqrt(-1)*(sgnom);
d=vrep(d,size(x,2),2);
y=ifft(X.*d);

%Take real part if real vector was input
if bool
  y=real(y);
end

function[]=anatrans_test

load solomon 
use solomon

x=[x y z];
z=anatrans(x,'periodic');

res=frac(vsum(abs(real(z)-x).^2,1),vsum(abs(z).^2,1));
reporttest('ANATRANS departure of real part from real-valued original less than 1/1000, Solomon Islands',allall(res<1/1000))

[zp,zn]=anatrans(z,conj(z),'periodic');
res=frac(vsum(abs(z-zp).^2,1),vsum(abs(z).^2,1));
reporttest('ANATRANS positive rotary part recovers input analytic signal, Solomon Islands',allall(res<1/1000))
res=vsum(abs(zn).^2,1);
reporttest('ANATRANS negative rotary part negligible for input analytic signal, Solomon Islands',allall(abs(zn)<1e-8))

z2=anatrans(x',2,'periodic');
reporttest('ANATRANS transpose orientation, periodic conditions',aresame(conj(z2'),z,1e-10))

z=anatrans(x,'periodic');
z2=anatrans(x',2,'periodic');
reporttest('ANATRANS transpose orientation, mirror conditions',aresame(conj(z2'),z,1e-10))

[z2,zp2,zn2,z4]=anatrans(x,z,conj(z),x,'periodic');
bool=aresame(z2,z)&&aresame(zp,zp2)&&aresame(zn,zn2)&&aresame(z4,z);
reporttest('ANATRANS multiple input arguments',bool)


