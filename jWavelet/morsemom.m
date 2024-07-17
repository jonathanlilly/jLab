function[m,n,k,l]=morsemom(p,ga,be)
%MORSEMOM  Frequency-domain moments of generalized Morse wavelets.
%
%   MORSEMOM is a low-level function called by several other Morse wavelet
%   functions.
%  
%   [MP,NP]=MORSEMOM(P,GAMMA,BETA) computes the Pth order frequency-
%   domain moment M and energy moment N of the lower-order generalized 
%   Morse wavelet specified by parameters GAMMA and BETA.
%
%   The Pth moment and energy moment are defined as
%
%           mp = 1/(2 pi) int omega^p  psi(omega)     d omega 
%           np = 1/(2 pi) int omega^p |psi(omega)|.^2 d omega 
%
%   respectively, where omega is the radian frequency.  These are evaluated
%   using the 'bandpass' normalization, which has max(abs(psi(omega)))=2.
%
%   The input parameters must either be matrices of the same size, or
%   some may be matrices and the others scalars.   
%
%   [MP,NP,KP,LP]=MORSEMOM(...) also returns the Pth order cumulant KP and
%   the Pth order energy cumulant LP.
%
%   Note that for very large BETA, the standard form of the moments fail
%   for numerical reasons, and one must use asymptotic forms.  These are
%   used by default when the argument to the gamma function exceeds 100 
%   and when then standard expressions yield non-finite or zero values. 
%
%   For details see
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   See also MORSEWAVE, MORSEDERIV, MOM2CUM.
%
%   'morsemom --t' runs some tests.
%
%   Usage:  mp=morsemom(p,ga,be);
%           [mp,np]=morsemom(p,ga,be);
%           [mp,np,kp,lp]=morsemom(p,ga,be);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2023 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(p, '--t')
    morsemom_test,return
end


m=morsemom1(p,ga,be);

if nargout>1
    n=frac(morseafun(ga,be).^2,morseafun(ga,2*be)).*frac(1,2.^(frac(2*be+1+p,ga))).*morsemom1(p,ga,2*be);
    n=frac(1,2.^(frac(1+p,ga))).*morsemom1(p,ga,2*be);
end

if nargout>2
     for i=0:maxmax(p)
         mcell{i+1}=morsemom1(i,ga,be);
     end
     kcell=mom2cum(mcell);
     if length(p)==1 
        k=kcell{p+1};
     else
        k=zeros(size(m));
        for i=1:length(m)
            k(i)=kcell{i};
        end
     end
end


if nargout>3
     for i=0:maxmax(p)
          ncell{i+1}=frac(1,2.^frac(1+i,ga)).*morsemom1(i,ga,2*be);
     end
     lcell=mom2cum(ncell);
     if length(p)==1 
        l=lcell{p+1};
     else
        l=zeros(size(m));
        for i=1:length(m)
            l(i)=lcell{i};
        end
     end
end



function[m]=morsemom1(p,ga,be)
m=morsemom1_standard(p,ga,be);

%correction for asymptotics 
if anyany(~isfinite(m))
    garg=frac(be+p+1,ga);
    if length(be)==1&&length(m)~=1
        be=be+zeros(size(m));
    end
    if length(p)==1&&length(m)~=1
        p=p+zeros(size(m));
    end
    if length(ga)==1&&length(m)~=1
        ga=ga+zeros(size(m));
    end
    bool=(~isfinite(m)|m==0)&(garg>100);
    %vsize(p,ga,be)
    m(bool)=morsemom1_asymptotic(p(bool),ga(bool),be(bool));
end

%morsef(ga,be+p)
%function[m]=morsemom1(p,ga,be)
%lnm=log(morseafun(ga,be))+log(frac(1,2*pi*ga))+gammaln(frac(be+p+1,ga));
%m=exp(lnm);

function[m]=morsemom1_asymptotic(p,ga,be)
fact=sqrt(frac(2,pi.*ga.*(be+p+1)));
m=fact.*frac(be+p+1,ga).^((p+1)./ga);

function[m]=morsemom1_standard(p,ga,be)
m=morseafun(ga,be).*morsef(ga,be+p);


function[f]=morsef(ga,be)
%MORSEF  Returns the generalized Morse wavelet first moment "f".
%
%   F=MORSEF(GAMMA,BETA) returns the normalized first frequency-
%   domain moment "F_{BETA,GAMMA}" of the lower-order generalized 
%   Morse wavelet specified by parameters GAMMA and BETA.

f=frac(1,2*pi*ga).*gamma(frac(be+1,ga));


function[]=morsemom_test

morsemom_test_expression;
morsemom_test_numerical;
morsemom_test_cumulant;
morsemom_test_energycumulant;
morsemom_test_alternate;
morsemom_test_ratio;

%function[]=morsemom_test_asymptotic
%A little plot to check asymptotic forms
%be=[4:1000]';
%m1=morsemom1_standard(-4,3,be);
%m2=morsemom1_asymptotic(-4,3,be);
%plot(be,[m1 m2]),ylog

%plot(morsemom(-4,3,be))
%ylog

function[]=morsemom_test_energycumulant 
ga=[3 0.5 9];be=[2.0000 15.1200 0.751];

[a,sigt,sigo]=morsebox(ga,be);
om=morsefreq(ga,be);
[m2,n2,k2,l2]=morsemom(2,ga,be);
reporttest('MORSEMOM second energy cumulant matches MORSEBOX',aresame(sqrt(l2./om./om),sigo,1e-6))

function[]=morsemom_test_expression
ga=(2:.1:10);
be=(1:.1:10);
[ga,be]=meshgrid(ga,be);

clear bool
for p=1:10
    [mp,np]=morsemom(p,ga,be);
    mp2=frac(exp(1)*ga,be).^(be./ga).*frac(1,2*pi*ga).*gamma(frac(be+1+p,ga));
    bool(p)=aresame(mp,mp2,1e-10);
end

reporttest('MORSEMOM analytic expression',all(bool))

function[]=morsemom_test_numerical

ga1=(2:2:12);
be1=(1:2:10);

%n=0;
omi=(0:.01:20)';
dom=omi(2)-omi(1);
                
[m1c,n1c,m2c,n2c]=vzeros(length(ga1),length(be1));
clear m1c n1c m2c n2c
for i=1:length(ga1)
    for j=1:length(be1)
        ompsi=morsefreq(ga1(i),be1(j));
        
        psi=morseafun(ga1(i),be1(j)).*(omi.*ompsi).^be1(j).*exp(-(omi.*ompsi).^ga1(i));
        %n=n+1;
        %psiall(:,n)=psi;
        
        m1c(j,i)=frac(dom,2*pi).*ompsi.^2.*vsum(omi.*psi,1);
        n1c(j,i)=frac(dom,2*pi).*ompsi.^2.*vsum(omi.*psi.^2,1);
        m2c(j,i)=frac(dom,2*pi).*ompsi.^3.*vsum(omi.^2.*psi,1);
        n2c(j,i)=frac(dom,2*pi).*ompsi.^3.*vsum(omi.^2.*psi.^2,1);
    end
end

[ga,be]=meshgrid(ga1,be1);
[m1,n1]=morsemom(1,ga,be);
[m2,n2]=morsemom(2,ga,be);

tol=2*10.^(-2);
bool=aresame(m1,m1c,tol)&&aresame(n1,n1c,tol)&&aresame(m2,m2c,tol)&&aresame(n2,n2c,tol);

reporttest('MORSEMOM numerical calculation',all(bool))

function[]=morsemom_test_cumulant

ga1=(2:2:12);
be1=(1:2:10);
[ga,be]=meshgrid(ga1,be1);

[mo,no,ko]=morsemom(0,ga,be);
[m1,n1,k1]=morsemom(1,ga,be);
[m2,n2,k2]=morsemom(2,ga,be);
[m3,n3,k3]=morsemom(3,ga,be);
    
kob=log(mo);
k1b=frac(m1,mo);
k2b=frac(m2,mo)-frac(m1,mo).^2;
k3b=2*frac(m1.^3,mo.^3)-3.*frac(m1,mo).*frac(m2,mo)+frac(m3,mo);


reporttest('MORSEMOM first three cumulants',aresame(ko,kob,1e-10)&&aresame(k1,k1b,1e-10)&&aresame(k2,k2b,1e-10)&&aresame(k3,k3b,1e-10))

function[]=morsemom_test_alternate
%Test alternate expressions for the energy moments

ga1=(2:2:12);
be1=(1:2:10);
[ga,be]=meshgrid(ga1,be1);

[mo,no]=morsemom(0,ga,be);
[m1,n1]=morsemom(1,ga,be);
[m2,n2]=morsemom(2,ga,be);
[m3,n3]=morsemom(3,ga,be);

nall(:,:,1)=no;
nall(:,:,2)=n1;
nall(:,:,3)=n2;
nall(:,:,4)=n3;

a=morseafun(ga,be);
nallb(:,:,1)=frac(a.^2,2.^((2*be+1)./ga)).*frac(1,2*pi*ga).*gamma(frac(2*be+1,ga));
nallb(:,:,2)=frac(a.^2,2.^((2*be+1+1)./ga)).*frac(1,2*pi*ga).*gamma(frac(2*be+1+1,ga));
nallb(:,:,3)=frac(a.^2,2.^((2*be+1+2)./ga)).*frac(1,2*pi*ga).*gamma(frac(2*be+1+2,ga));
nallb(:,:,4)=frac(a.^2,2.^((2*be+1+3)./ga)).*frac(1,2*pi*ga).*gamma(frac(2*be+1+3,ga));

reporttest('MORSEMOM alternate expressions for energy moments',aresame(nall, nallb, 1e-6));

[m3,n3,k3,l3]=morsemom(3,ga,be);
l3b=frac(n3,no)-3*frac(n1.*n2,no.^2)+2*frac(n1.^3,no.^3);

reporttest('MORSEMOM alternate expression for third energy cumulant',aresame(l3,l3b,1e-6));

function[]=morsemom_test_ratio

ga1=(2:.5:12);
be1=(1:.5:10);
[ga,be]=meshgrid(ga1,be1);


[mo,no]=morsemom(0,ga,be);
[m1,n1]=morsemom(1,ga,be);
[m2,n2]=morsemom(2,ga,be);

%momrat=frac(gamma(frac(2*be+3,ga)).*gamma(frac(2*be+1,ga)),squared(gamma(frac(2*be+2,ga))))-1;
momrat=exp(gammaln(frac(2*be+3,ga))+gammaln(frac(2*be+1,ga))-2*gammaln(frac(2*be+2,ga)))-1;
momrato=(n2./no)./squared(n1./no)-1;

reporttest('MORSEMOM moment ratio expression matches calculation',aresame(momrat,momrato, 1e-6));

% ga1=(1:.1:12);
% be1=(1:.1:10);
% [ga,be]=meshgrid(ga1,be1);
% [m3,n3,k3,l3]=morsemom(3,ga,be);
% 
% 
% ga1=(3:.001:4);
% be1=(1:.1:40);
% [ga,be]=meshgrid(ga1,be1);
% [m3,n3,k3,l3]=morsemom(3,ga,be);
% [m,index]=min(abs(l3),[],2);
% figure,plot(be1,ga1(index))

