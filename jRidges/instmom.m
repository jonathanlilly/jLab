function[varargout]=instmom(varargin)
%INSTMOM  Univariate and multivariate instantaneous moments.
%
%   [A,OMEGA,UPSILON]=INSTMOM(X), where X is an analytic signal, computes 
%   the amplitude A, instantaneous *radian* frequency OMEGA, and 
%   instantaneous bandwidth assuming a unit sample rate. 
%
%   X is an array with the first dimension being "time".  Thus, X can be a 
%   matrix of analytic signals oriented as column vectors, or a 2- or 3-D 
%   wavelet transform such as output by WAVETRANS.
%
%   The output arrays are the same size as X. 
%
%   The instantaneous frequency, bandwidth, and curvature are defined as
%
%            A    = abs X
%          OMEGA  = d/dt Im ln X = d/dt arg X
%         UPSILON = d/dt Re ln X = d/dt ln abs X  
%  
%   where i=SQRT(-1) as usual.
%
%   INSTMOM(X,DIM) computes the moments along dimension DIM, instead of 
%   the default of computing the moments along the rows (DIM=1).
%
%   For details, see 
%  
%       Lilly & Olhede (2010), "Bivariate instantaneous frequency and 
%            bandwidth", IEEE Trans. Sig. Proc., 58 (2), 591--603.
%   _____________________________________________________________________
%   
%   Sample interval
%
%   INSTMOM(DT,...) uses sample interval DT, where DT is a scalar, for 
%   computing time derivatives.  DT=1 is the default.
%   _____________________________________________________________________
% 
%   Joint instantaneous moments
%
%   INSTMOM can also calculate the joint instananeous moments of 
%   multivariate signals, as defined in Lilly and Olhede (2010).
%
%   [JA,JOMEGA,JUPSILON]=INSTMOM(X1,X2,...,XN,DIM) returns the *joint*
%   instantaneous moments calculated across the N signals X1,X2,... XN, 
%   based on the univariate instantaneous moments along dimension DIM.
%
%   The joint instantaneous amplitude JA is the root-mean-square of the 
%   component amplitudes across dimension JDIM, while JOMEGA is power-
%   weighted average of the component instantaneous frequencies.
%
%   For details and for the definition of the joint instantaneous bandwidth 
%   JUPSILON, see Lilly and Olhede (2010).
%
%   [JA,JOMEGA,JUPSILON]=INSTMOM(X,DIM,JDIM) also works, where the joint
%   instantaneous moments are calculated across dimensions JDIM of X.
%
%   The joint instantaneous moments JA, JOMEGA, and JUPSILON then have the 
%   same size as X, except along dimension JDIM where they have only one 
%   entry.  Note that DIM is no longer optional when JDIM is used.
%   _____________________________________________________________________
% 
%   Instantaneous curvature
%
%   INSTMOM can also return the next-higher order instantaneous moment, 
%   which is more rarely encountered. 
%
%   [A,OMEGA,UPSILON,XI]=INSTMOM(X) returns the instantaneous curvature XI,
%   defined as
%
%            XI   = d^2/dt^2 abs X / abs X + i d^2/dt^2 arg X
%                 = UPSILON^2 + d/dt UPSILON + i d/dt OMEGA
%
%   Similarly [JA,JOMEGA,JUPSILON,JXI]=INSTMOM(X,DIM,JDIM) returns the
%   joint instantaneous curvature JXI for the multivariate signal X.
%    
%   For details on the univariate and joint instantaneous curvature, see 
%
%       Lilly and Olhede (2012a), "Analysis of modulated multivariate 
%            oscillations", IEEE Trans. Sig. Proc., 60 (2), 600--612. 
%   _____________________________________________________________________
%
%   Boundary conditions
%
%   The first and last points must be treated differently, as the central 
%   difference is not defined there.  Three different methods can be used.
%
%   INSTMOM(...,STR) specifies the method: STR= 'endpoint' (the default),
%   'periodic', or 'nans'.  See VDIFF for details.   
%   _____________________________________________________________________
%
%   'instmom --f' generates some sample figures.
%
%   Usage: [a,om]=instmom(x);
%          [a,om,up,xi]=instmom(dt,x);
%          [a,om,up,xi]=instmom(dt,x,dim);
%          [a,om,up,xi]=instmom(dt,x1,x2,x3,x4,dim);
%          [a,om,up,xi]=instmom(dt,x,dim,jdim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2018 J.M. Lilly --- type 'help jlab_license' for details


%   _____________________________________________________________________
%
%   Higher-order modulation functions
%
%   [OMEGA,RHO1,...RHON]=INSTMOM(X) also outputs the higher-order 
%   instantaneous modulation functions.  RHO1 is identical to the
%   bandwidth, and RHO2 is called the curvature.
%
%   For details see 
%
%      Lilly and Olhede (2010).  On the analytic wavelet transform.
%
%   Note that the modulation functions are defined in their non-normalized 
%   form, that is, not divided by powers of the instananeous frequency, 
%   unlike in Lilly and Olhede (2010).
%
%   _____________________________________________________________________


if strcmpi(varargin{1}, '--t')
    instmom_test,return
elseif strcmpi(varargin{1}, '--f')
    type makefigs_instmom
    makefigs_instmom;
    return
end

if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else
    dt=1;
end

        
str='endpoint';
dim=1;
jdim=[];
for i=1:4
    if ischar(varargin{end})
        str=varargin{end};
        varargin=varargin(1:end-1);
    else
        if length(varargin{end})==1
            if length(varargin{end-1})==1    
                dim=varargin{end-1};
                jdim=varargin{end};
                varargin=varargin(1:end-2);
            else
                dim=varargin{end};
                varargin=varargin(1:end-1);
            end
        end
    end
end




N=length(varargin);
if N==1
    x=varargin{1};
else
    jdim=lnsd(varargin{1})+1;
    sizex1=size(varargin{1});
    sizex1=sizex1(1:jdim-1);
    x=zeros([sizex1 N]);
    for i=1:N
         x=vindexinto(x,varargin{i},i,jdim);
    end
end

%Note: there are different ways to reasonably define the instantaneous
%frequency as a vdiff of the signal.  These differ greatly actually. 
%Diffing the unwrapped angle seems to give the best results.  If you take 
%the Im of the diffed log, you will break the wavelet phase algorithm!

%You definitely don't want to just do diff(x).  This differences the real
%and imaginary parts separately, which are themselves rapidly varying.
  
om=vdiff(unwrap(angle(x),[],dim),dim,str)./dt;
up=vdiff(log(abs(x)),dim,str)./dt; 
eta=om-sqrt(-1)*up;

varargout{1}=abs(x);
varargout{2}=om;
nmax=nargout-2;

%I prefer this version, though difference is minor
%if nmax>=1
%   varargout(3)=frac(1,abs(x)).*vdiff(vdiff(abs(x),1,str),1,str)./dt./dt+sqrt(-1)*vdiff(om,1,str)./dt; 
%end

if nmax>=1
    etadiff=cell(nmax,1);   
    polyargs=cell(nmax,1);
    bcell=cell(nmax,1);

    etadiff{1}=vdiff(eta,dim,str)./dt;
    polyargs{1}=up;

    for n=2:nmax
        etadiff{n}=vdiff(etadiff{n-1},dim,str)./dt;
        polyargs{n}=sqrt(-1)*etadiff{n-1};
    end

    temp=bellpoly(polyargs);
    for n=1:length(temp)
        varargout{n+2}=temp{n};
    end
end

%When JDIM is input, this tells me to output joint quantities
if ~isempty(jdim)
    if nargout==1
        varargout{1}=jointmom(x,jdim);    
    elseif nargout==2
        [varargout{1},varargout{2}]=jointmom(x,varargout{2},jdim);
    elseif nargout==3
        [varargout{1},varargout{2},varargout{3}]=jointmom(x,varargout{2},varargout{3},jdim);
    elseif nargout==4
        [varargout{1},varargout{2},varargout{3},varargout{4}]=jointmom(x,varargout{2},varargout{3},varargout{4},jdim);
    end
end


function[varargout]=jointmom(varargin)
%JOINTMOM  Joint instantaneous frequency, bandwidth, and curvature.
%
%   OMEGAX=JOINTMOM(X,OMEGA,DIM) where X is an array of multiple analytic
%   signals and OMEGA is an array of their instantaneous frequencies, gives
%   the joint instaneous frequency OMEGAX averaged over dimension DIM.
%
%   X is presumed to have time oriented in rows.  The output matrix OMEGAX
%   are the same size as X except along dimension DIM, where the output
%   matrices will have length one.  X and OMEGA are the same size.
%   
%   [OMEGAX,UPSILONX]=JOINTMOM(X,UPSILON,DIM) also works, where UPSILON 
%   are the individual bandwidths, and UPSILONX is the joint quantity. 
%   UPSILON is the same size as X and OMEGA.
%
%   Finally [OMEGAX,UPSILONX,XIX]=JOINTMOM(X,UPSILON,XI,DIM) also returns 
%   the joint instantaneous curvature XIX given individual curvatures XI.
%
%   For details, see 
%  
%       Lilly & Olhede (2010), "Bivariate instantaneous frequency and 
%            bandwidth", IEEE Trans. Sig. Proc., 58 (2), 591--603.
% 
%   See also INSTMOM.
%
%   Usage: [a,om]=jointmom(x,om,dim);
%          [a,om,up]=jointmom(x,om,up,dim);
%          [a,om,up,xi]=jointmom(x,om,xi,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2014 J.M. Lilly --- type 'help jlab_license' for details
 


dim=varargin{end};
x=varargin{1};
varargin=varargin(2:end-1);


if ~isempty(varargin)
    om=varargin{1};
end
if length(varargin)>1
    upsilon=varargin{2};
end
if length(varargin)>2
   xi=varargin{3};
end

varargout{1}=sqrt(vmean(abs(x).^2,dim));
if nargout>1
    varargout{2}=vmean(om,dim,squared(x));
    ombar=vrep(varargout{2},size(x,dim),dim);
end
if nargout>2
    varargout{3}=sqrt(vmean(abs(upsilon+sqrt(-1)*(om-ombar)).^2,dim,squared(x)));
end
if nargout>3
    varargout{4}=sqrt(vmean(abs(xi+2*sqrt(-1)*upsilon.*(om-ombar)-(om-ombar).^2).^2,dim,squared(x)));
end
if nargout>4
     error('Sorry, INSTMOM only outputs the first two deviation vectors for joint moments.')
end


function[]=instmom_test

load npg2006
use npg2006

[a1,om1,up1,c1]=instmom(cv);
[a2,om2,up2,c2]=instmom(permute(cv,[2 1]),2);

reporttest('INSTMOM moments computed along different dimensions match', aresame([om1 up1 c1],[om2(:) up2(:) c2(:)])); 
jointmom_test;

function[]=jointmom_test
load solomon
use solomon

[x,z]=anatrans(x,z);

[ax,omx,upx]=instmom(x);
[az,omz,upz]=instmom(z);

ombar=frac(abs(x).^2.*omx+abs(z).^2.*omz,abs(x).^2+abs(z).^2);
upbar=sqrt(frac(abs(x).^2.*(upx.^2+(omx-ombar).^2)+abs(z).^2.*(upz.^2+(omz-ombar).^2),abs(x).^2+abs(z).^2));

[a,om,up]=instmom([x z],1,2);

reporttest('INSTMOM joint moments using Solomon Islands frequency',aresame(om,ombar,1e-8))
reporttest('INSTMOM joint moments using Solomon Islands bandwidth',aresame(up,upbar,1e-8))

[a2,om2,up2]=instmom(x,z,1);

reporttest('INSTMOM joint moments alternate form',aresame([a om up],[a2 om2 up2]))



