function[t,om]=morseregion(varargin)
%MORSEREGION  Generalized Morse wavelet time-frequency concentration region.
%
%   MORSEREGION computes two different type of time-frequency concentration
%   regions for the generalized Morse wavelets.
%
%   This function is part of 'element analysis' described in Lilly (2017), 
%   "Element analysis: a wavelet-based method for analyzing time-localized
%   events in noisy time series", submitted.  Available at www.jmlilly.net.
%   __________________________________________________________________
%
%   Localization regions
%
%   [T,F]=MORSEREGION(A,GAMMA,BETA,FO) returns the localization region of 
%   the generalized Morse wavelets specified by GAMMA and BETA centered on
%   *radian* frequency FO.  The nonnegative number A sets the region area.
%
%   T and F are then a parametric curve such that PLOT(T,F) shows the shape
%   of the time-frequency localization region. 
%
%   The localization regions are based on a reconstruction of the signal 
%   from a limited inversion region, see Daubechies and Paul (1988) and 
%   Olhede and Walden (2002, 2003a) for details.  
% 
%   Note that the definition of area used here is one-half of that used in
%   Olhede and Walden (2002, 2003a), as we only count positive frequencies.
%   __________________________________________________________________
%
%   Regions of influence
%
%   [T,F]=MORSEREGION(LAMBDA,GAMMA,BETA,MU,FRHO) with five input arguments 
%   returns the contour at which the energy-normalized wavelet transform of
%   another wavelet falls to a fraction LAMBDA of its peak value.
%
%   The wavelet we are taking the transform of has parameters GAMMA and MU,
%   and is characterized by a scale frequency FRHO.  We then take the 
%   wavelet transform of this wavelet with a (GAMMA,BETA) wavelet.  Note
%   that this region assumes a 1/S and not a 1/SQRT(S) normalization.
%
%   T and F are then a parametric curve such approximating the region where
%   this transform takes on a fraction LAMBDA of its peak value.  This 
%   approximation is formed analytically through a Taylor series expansion.
%
%   Note that FRHO is the frequency of the analyzed wavelet, not the 
%   location of the maximum within the transform, which is given by
%
%       FO=FRHO*MORSEFREQ(GAMMA,BETA)/MORSEFREQ(GAMMA,MU)*
%                                      ((MU+1)./BETA)^(1/GAMMA).
%
%   For futher details, see Lilly (2017).
%   __________________________________________________________________
%
%   'morseregion --f' generates a sample figure.
%
%   Usage: [t,f]=morseregion(A,ga,be,fo);
%          [t,f]=morseregion(lambda,ga,be,mu,frho);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2017 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--f')
   type makefigs_morseregion
   makefigs_morseregion;
   return
elseif strcmpi(varargin{1},'--t')
    morsecfun_test;return
end

str='band';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

A=varargin{1};
ga=varargin{2};
be=varargin{3};
if length(varargin)==5
    mu=varargin{4};
else
    mu=[];
end
frho=varargin{end};

%if nargin==3
%    frho=(fm./fe);
%end
%These are the same, apparently
%t=frac(c2,a.^(1-1./ga)).*real(sqrt(2*a.*C-1-a.^2));

N=1000;
if isempty(mu)
    arrayify(A,ga,be,frho);
    [t,om]=morseregion_localization(N,A,ga,be,frho);
else
    arrayify(A,ga,be,mu,frho);
    [t,om]=morseregion_influence(N,A,ga,be,mu,frho,str);
end

function[t,om]=morseregion_influence(N,C,ga,be,mu,frho,str)

[fm,fe]=morsefreq(ga,be);
[t,s]=vzeros(N,length(C),'nan');

for i=1:length(C)
    [mp,np,k2,lp]=morsemom(2,ga(i),be(i)+mu(i));
    
    if strcmpi(str(1:3),'ban')
        smax=frac(be(i),mu(i)+1).^(1./ga(i));
        B=frac(smax.^(be(i)),(smax.^ga(i)+1).^((be(i)+mu(i)+1)./ga(i)));
        
        slow=(C.*B).^(1./be(i));
        shigh=frac(1,C.*B).^(1./(mu(i)+1));
        so=logspace(log10(slow),log10(shigh),N)';
        x1=be(i).*log(so);
    elseif strcmpi(str(1:3),'ene')
        smax=frac(be(i)+1/2,mu(i)+1/2).^(1./ga(i));
        B=frac(smax.^(be(i)+1/2),(smax.^ga(i)+1).^((be(i)+mu(i)+1)./ga(i)));
        
        slow=(C(i).*B).^(1./(be(i)+1/2));
        shigh=frac(1,C(i).*B).^(1./(mu(i)+1/2));
        so=logspace(log10(slow),log10(shigh),N)';
        x1=(be(i)+1/2).*log(so);
    end

    %First I compute t and omega for the beta and gamma 
    fact=sqrt(2)*frac((so.^ga(i)+1).^(1./ga(i)),sqrt(k2));
    x2=frac(be(i)+mu(i)+1,ga(i)).*log(so.^ga(i)+1);
    ti=(fact.*sqrt(-log(C(i))-log(B)+x1-x2));
    si=so;
    
    %Remove points for which the solution is not real-valued
    vindex(ti,si,find(abs(imag(ti))==0),1);
    if ~isempty(ti)
        %Flip-concatenate and add the last point to the beginning
        si=[si;flipud(si)];
        ti=[-ti;flipud(ti)];
        
        ti=[ti(end);ti];
        si=[si(end);si];
        
        %Interpolate to length N
        t(:,i)=interp1([0:length(ti)-1]'./(length(ti)-1),ti,[0:N-1]'./(N-1),'linear');
        s(:,i)=interp1([0:length(si)-1]'./(length(si)-1),si,[0:N-1]'./(N-1),'linear');
    end
end
rho=frac(morsefreq(ga,mu),frho);
t=t.*rho;
om=frac(fm,s*rho);

function[t,om]=morseregion_localization(N,A,ga,be,fs)

[fm,fe]=morsefreq(ga,be);
C=morsecfun(A,ga,be);

r=((2*be)+1)./ga;
c1=2.^(-1./ga).*frac(gamma(r+1./ga),gamma(r));
c2=frac(be,ga).*2.^(1./ga).*frac(gamma(r-1./ga),gamma(r));

coeff=gamma(r+(1./ga)).*2.^(-1./ga);
ommin=coeff./(gamma(r).*(C+sqrt(C.^2-1)).^(1./ga));
ommax=coeff./(gamma(r).*(C-sqrt(C.^2-1)).^(1./ga));

om=zeros(N/2,length(C));
for i=1:length(C)
    %om(:,i)=linspace(ommin(i),ommax(i),N/2)';
    om(:,i)=logspace(log10(ommin(i)),log10(ommax(i)),N/2)';
end

vtranspose(ga,be,c1,c2,C);
vrep(ga,be,c1,c2,C,N/2,1);
a=frac(c1,om).^ga;
b=frac(a.^(1-1./ga),c2);
t=real(frac(sqrt(2*a.*C-1-a.^2),b));

om=[om;flipud(om)];
t=[-t;flipud(t)];

vtranspose(fs,fm);
vrep(fs,fm,N,1);
om=om.*(fs./fm);
t=t./(fs./fm);

function[c]=morsecfun(a,ga,be)
%MORSECFUN  Morse wavelet "C"-function.
%
%   C=MORSECFUN(A,GAMMA,BETA) returns the value of the generalized Morse
%   wavelet "C"-parameter in terms of the area of concentration A and the
%   GAMMA and BETA parameters. 
%
%   The input parameters may either be arrays of the same size, or some
%   may be arrays and the others scalars.  
%
%   MORSECFUN uses the formula of Olhede and Walden (2002), "Generalized 
%   Morse Wavelets", for the area of "D_{C,BE,GA}" given at the bottom 
%   right of page 2664.  
%
%   Note the area used by MORSECFUN and MORSEAREA is a "one-sided" version,
%   differing by Olhede and Walden's by a factor of 1/2.  
%
%   See also MORSEAREA, MORSEREGION.
%
%   'morsecfun --t' runs a test.
% 
%   Usage: C = morsecfun(A,ga,be);  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2017 J.M. Lilly --- type 'help jlab_license' for details        
  

arrayify(a,ga,be);

r=((2*be)+1)./ga;
fact=2.*pi.*gamma(r+1-(1./ga)).*gamma(r+(1./ga))./(ga.*gamma(r).^2);
c=2*a./fact+1;
%Note this differs from Olhede and Walden by a factor of 2, since
%I define area to be different by a factor of 1/2

function[b]=morsecfun_test
be=1;
ga=1;
A=[10 150]';

C= morsecfun(A,ga,be);
A2=morsearea(C,ga,be);

tol=1e-10;
b=aresame(A,A2,tol);
reporttest('MORSECFUN inverts MORSEAREA',b);


function[]=morseregion_fig
 
ga=2;
be=2;
mu=2;
alpha=0.9;
smax=frac(be,mu+1).^(1./ga);
B=frac(smax.^(be+1/2),(smax.^ga+1).^((be+mu+1)./ga));
slow=(alpha.*B).^(1./(be+1/2));
shigh=frac(1,alpha.*B).^(1./(mu+1/2));
so=[0.0001:0.0001:100]';
y=[so.^(be+1/2) alpha.*B.*(1+so.^ga).^((be+mu+1)./ga)];
figure,plot(log10(so),log10(y)),vlines(log10(slow)),vlines(log10(shigh))
%This illustrates the bounds used for the region of influence calculation



figure
ga1=[ 1 3 9 27];
be1=[1 3 9 27];
[ga,be]=meshgrid(ga1,be1);
[t,f]=morseregion(2,ga,be,1);
vcolon(ga,be);

subplot(211)
h=plot(t.*oprod(ones(size(t(:,1))),1./(be(:))),f,'k');
linestyle(h(be==2),'k-.')
linestyle(h(be==4),'k-')
linestyle(h(be==8),'k')
title('Morse wavelet localization regions')
xlabel('Normalized Time')
ylabel('Normalized Frequency, Linear Axis'),
linestyle default
hlines(1,':'),vlines(0,':')

subplot(212)
h=plot(t.*oprod(ones(size(t(:,1))),1./(be(:))),f,'k');
xlabel('Normalized Time')
ylabel('Normalized Frequency, Logarithmic Axis'),
linestyle default
ylog
hlines(1,':'),vlines(0,':')

