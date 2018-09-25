function[nui,om2i]=morlfreq(omi)
%MORLFREQ  Compute Morlet wavelet carrier frequency given peak frequency.
%
%   FNU=MORLFREQ(FMAX), where FMAX is the desired peak frequency of a 
%   Morlet wavelet, returns the appropriate carrier wave frequency FNU.
%
%   The Morlet wavelet carrier wave and peak frequency differ on account
%   of the correction term which ensures zero mean.
%
%   Note that FMAX and FNU are both *radian* as in cos(omega t) and not
%   cyclic as in in cos(2 pi f t).
%
%   [FNU,FMIN]=MORLFREQ(FMAX) also returns the (radian) frequency at which
%   the Morlet wavelet obtains a minimum, which occurs for FMIN<0.
%
%   For details see Appendix A of
%   
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%          wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   'morlfreq --t' runs a test.
%   'morlfreq --f' generates a sample figure.
%
%   Usage: fnu=morlfreq(fmax);
%          [fnu,fmin]=morlfreq(fmax);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2013 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(omi, '--t')
    morlfreq_test,return
end
if strcmpi(omi, '--f')
    type makefigs_morlfreq
    makefigs_morlfreq;
    return
end
 
sizeom=size(omi);
omi=omi(:);

N=10000;

%This gives a tighter spacing near one than near zero
nutilde=log(log(linspace(exp(exp(1./1e9)),exp(exp(1-1./1e9)),N)));
%nutilde=linspace(1./1e9,1-1./1e9,N);

om=(sqrt(-frac((log(1-nutilde)),nutilde)));
nu=nutilde.*om;

bool=(omi>=maxmax(om));
index1=find(bool);
if ~isempty(index1)
    nui(index1)=omi(index1);
end
index2=find(~bool);
if ~isempty(index2)
    nui(index2)=interp1(om,nu,omi(index2),'linear');
end

%I chose this grid via trial and error
nutilde=1./[linspace(1./10,10,N) logspace(log10(10+0.0000001),log10(1000),N)];
om2=-(sqrt(frac((log(1+nutilde)),nutilde)));
nu2=(nutilde).*abs(om2);

om2i=zeros(size(nui));%figure,plot(nui),maxmax(nu2)
bool=(nui>=maxmax(nu2));
index1=find(bool);
if ~isempty(index1)
    om2i(index1)=nan;
end
index2=find(~bool);%length(index2)
if ~isempty(index2)
    om2i(index2)=interp1(nu2,om2,nui(index2),'linear');
end
%x=om2i-nui-om2i.*exp(-om2i.*nui);
%figure,plot(x)

nui=reshape(nui,sizeom);
om2i=reshape(om2i,sizeom);

function[]=morlfreq_test
om=linspace(0,4,1000);
[nu,fneg]=morlfreq(om);

psi1a=exp(-frac(1,2).*(nu-om).^2).*(-(om-nu))+exp(-frac(1,2).*(nu.^2+om.^2)).*om;
reporttest('MORLFREQ first derivative vanishes at peak',aresame(psi1a,0*psi1a,1e-4))

omneg=fneg*2*pi;
psi1b=exp(-frac(1,2).*(nu-omneg).^2).*(-(omneg-nu))+exp(-frac(1,2).*(nu.^2+omneg.^2)).*omneg;

reporttest('MORLFREQ first derivative vanishes at negative peak',aresame(psi1b,0*psi1b,1e-4))


