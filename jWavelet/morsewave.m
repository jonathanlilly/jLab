function[varargout]=morsewave(varargin)
% MORSEWAVE  Generalized Morse wavelets of Olhede and Walden (2002). 
%   _______________________________________________________________________
%
%   *|* morsewave.png --- Figure illustrating the generalized Morse 
%   wavelets. Type 'jhelp morsewave' to view this image. *|*
%   _______________________________________________________________________
%
%   PSI=MORSEWAVE(N,GAMMA,BETA,FS) returns an  N x LENGTH(FS) array PSI 
%   which contains time-domain versions of the generalized Morse wavelet
%   specified by GAMMA and BETA, concentrated at frequencies FS.
%  
%   The vector FS specifically denote the *radian* frequencies at which 
%   the Fourier transform of the wavelets reach their maximum amplitudes.  
%
%   A set of frequencies appropriate for analyzing a given length time 
%   series can be easily chosen using MORSESPACE.
%
%   Note that the wavelets are centered at the midpoint in time, that is, 
%   row number ROUND(SIZE(PSI,1)/2).  FS assumes a unit sample rate. 
%
%   [PSI,PSIF]=MORSEWAVE(...) optionally returns a frequency-domain version
%   PSIF of the wavelets.  PSIF is the same size as PSI.
%   _________________________________________________________________
%
%   Normalization
%
%   MORSEWAVE supports two kinds of normalization for the wavelets.
%
%   MORSEWAVE(...,'bandpass') uses "bandpass normalization", meaning that
%   the FFT of the wavelet has a peak value of 2 for all frequencies FS. 
%
%   MORSEWAVE(...,'energy') uses the unit energy normalization.  The time-
%   domain wavelet energy SUM(ABS(PSI).^2,1) is then always unity. 
%
%   The bandpass normalization corresponds to having 1/S in the time-domain
%   wavelet transform defintion, where S is the scale, while the unit 
%   energy normalization corresponds to 1/SQRT(S).
%
%   MORSEWAVE uses bandpass normalization by default.
%   _________________________________________________________________
%   
%   Multiple orthogonal wavelets 
%
%   MORSEWAVE can compute multiple orthogonal versions of the generalized
%   Morse wavelets, characterized by the order K.
%
%   PSI=MORSEWAVE(N,K,GAMMA,BETA,FS) with a fifth numerical argument K
%   returns an N x LENGTH(FS) x K array PSI which contains time-domain 
%   versions of the first K orthogonal generalized Morse wavelets.
%
%   These K different orthogonal wavelets have been employed in 
%   multiwavelet polarization analysis, see Olhede and Walden (2003a,b).
%
%   Again either bandpass or energy normalization can be applied.  With
%   bandpass normalization, all wavelets are divided by a constant, setting
%   the peak value of the first frequency-domain wavelet equal to 2.
%   _________________________________________________________________
%
%   Background
%
%   For further details on generalized Morse wavelets, see the following 
%   publications. 
%
%     Lilly and Olhede (2012), Generalized Morse wavelets as a superfamily
%        of analytic wavelets. IEEE Trans. Sig. Proc., 60 (11), 6036--6041.
%
%     Lilly and Olhede (2009),  Higher-order properties of analytic 
%         wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%     Olhede and Walden (2002),  Generalized Morse Wavelets. IEEE Trans.
%         Sig. Proc., 50 (11), 2661--2670.
%   _________________________________________________________________
%
%   'morsewave --t' runs a test.
%   'morsewave --f' generates some sample figures.
%
%   Usage: psi=morsewave(N,ga,be,fs);
%          [psi,psif]=morsewave(N,ga,be,fs,'bandpass');
%          [psi,psif]=morsewave(N,K,ga,be,fs,'energy');
%          [psi,psif]=morsewave(N,K,ga,be,fs,'bandpass');
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2016 J.M. Lilly and F. Rekibi
%                         --- type 'help jlab_license' for details  
 

%   _________________________________________________________________
%
%   The zero beta case
%
%   It is unlikely that you will need to read this section, which is 
%   why the comment is hidden.  It describes a feature mostly used for 
%   testing purposes.
%
%   For BETA equal to zero, the generalized Morse wavelets describe
%   a non-zero-mean function which is not in fact a wavelet. 
%
%   Only 'bandpass' normalization is supported for this case.
%
%   In this case the frequency speficies the half-power point of the
%   analytic lowpass filter.  
%
%   The frequency-domain definition of MORSEWAVE is not necessarily 
%   a good way to compute the zero-beta functions, however.  You will
%   probably need to take a very small DT.
%   _________________________________________________________________


%
%   Primary and edge wavelets
%
%   PSI=MORSEWAVE(N,GAMMA,BETA,FS,'edge') returns an N x LENGTH(FS) x 2 
%   array PSI containing two different wavelets at each frequency. 
%
%   The first page of PSI, PSI(:,:,1), contains the standard generalized 
%   Morse wavelet, which we will call the 'primary' wavelet.
%
%   The second page of PSI, PSI(:,:,2), is called the 'edge' wavelet, and
%   is orthogonal to the primary wavelet.  The edge wavelet is formed by 
%   taking the first frequency-domain derivative of the primary wavelet.  
%
%   A transform with edge wavelet can be compared against a transform with
%   the primary in order to assess the statistical significance of wavelet
%   ridge analysis, see RIDGEWALK for details.
%
%   This option supports both 'bandpass' and 'energy' normalizations.
%
%   Under the 'energy' normalization, both the primary and edge wavelets
%   are set to have unit energy.  Under the 'bandpass' normalization, both
%   wavelets are subsequently divided by a constant such that the maximum 
%   value of the primary wavelet in the frequency domain is equal to 2.
%   _________________________________________________________________
%
%   The first page of PSI, PSI(:,:,1), is again the first or standard
%   generalized Morse wavelet.  Note however, that the K=2 wavelet is
%   different from the edge wavelet discussed above, because the procedure 
%   for orthogonalization in the two cases is somewhat different. 

%   This behavior is deprecated  
%
%   Sample rate
%
%   MORSEWAVE(DT,N,K,GAMMA,BETA,F) specifies the sample rate to be DT.
%   DT is optional with a default value of unity.
%   _________________________________________________________________

if strcmpi(varargin{1},'--t')
  morsewave_test;return
end
if strcmpi(varargin{1},'--f')
   type makefigs_morsewave
   makefigs_morsewave;
   return
end

%Contingency for WAVETRANS inputting cell array
if iscell(varargin{2})
    temp=varargin{2};
    for i=1:length(temp)
        varargin{1+i}=temp{i};
    end
end

str='ban';  %Bandpass default
fam='fir';  %First family default
dom='tim';  %Time computation default
K=1;        %Only the first wavelet by default
for i=1:3
    if ischar(varargin{end})||isempty(varargin{end})
        temp=lower(varargin{end});
        if length(temp)>=3
            temp=temp(1:3);
            if strcmpi(temp,'fir')||strcmpi(temp,'edg')
                fam=temp;
            elseif strcmpi(temp,'ban')||strcmpi(temp,'ene')
                str=temp;
            elseif strcmpi(temp,'tim')||strcmpi(temp,'fre')
                dom=temp;
            end
        end
        varargin=varargin(1:end-1);
    end
end

dt=1;
N=varargin{1};
if length(varargin)==5
    K=varargin{2};
end
ga=varargin{end-2};
be=varargin{end-1};
fs=varargin{end};

if strcmpi(fam,'edg')&&K==1
   K=2;
end

%ga,be,K,N
% if be==0&&strcmpi(str,'ene')
%     str='ban';
%     disp('For BETA=0, energy normalization is not defined.  Using bandpass normalization.')
% end

x=zeros(N,length(fs),K);
X=zeros(N,length(fs),K);

for n=1:length(fs)
    [X(:,n,:),x(:,n,:)]=morsewave1(N,K,ga,be,abs(fs(n)),str,fam);
    if fs(n)<0
        if ~isempty(x)
            x(:,n,:)=conj(x(:,n,:));
        end
        X(2:end,n,:)=flip(X(2:end,n,:),1);
    end
end

varargout{1}=x;
varargout{2}=X;

function[X,x]=morsewave1(N,K,ga,be,fs,str,fam)
  
x=[];
fo=morsefreq(ga,be);
fact=fs./fo;
om=2*pi*linspace(0,1-1./N,N)'./fact;

if strcmpi(str,'ene')
    if be==0
        psizero=exp(-om.^ga);
    else
        psizero=exp(be.*log(om)-om.^ga);
    end
else
    if be==0
        psizero=2*exp(-om.^ga);
    else
        %Alternate calculation to cancel things that blow up
        psizero=2*exp(-be.*log(fo)+fo.^ga+be.*log(om)-om.^ga);
    end
end
%figure,plot(psizero)
psizero(1)=1/2*psizero(1); %Due to unit step function
%Ensure nice lowpass filters for beta=0;
%Otherwise, doesn't matter since wavelets vanishes at zero frequency


vswap(psizero,nan,0);

if strcmpi(fam,'fir')
   X=morsewave_first_family(fact,N,K,ga,be,om,psizero,str); 
elseif strcmpi(fam,'edg')
   X= morsewave_second_family(fact,N,K,ga,be,om,psizero,str); 
end

X=vswap(X,inf,0);

ommat=vrep(vrep(om,size(X,3),3),size(X,2),2);
Xr=X.*(-1).^[0:N-1]'; %ensures wavelets are centered
if ~mod(N, 2)
  Xr(N/2+1, :)/=2; %ensures proper wavelet decay and analyticity
endif
%figure,plot(vrep(om,size(X,2),2).*(N+1)/2*fact)

x=ifft(Xr);
%Note to self, fft and ifft do not exactly invert one another

function[psif]=morsewave_first_family(fact,N,K,ga,be,om,psizero,str)

r=(2*be+1)./ga;
c=r-1;
L=0*om;
index=(1:floor(N/2)+1);
psif=zeros(length(psizero),1,K);

for k=0:K-1
  %Log of gamma function much better ... trick from Maltab's ``beta''
  if strcmpi(str,'ene')
        A=morseafun(k+1,ga,be,str);
        coeff = sqrt(1./fact)*A;
  elseif strcmpi(str,'ban')
        if be~=0
            coeff=sqrt(exp(gammaln(r)+gammaln(k+1)-gammaln(k+r)));
        else
            coeff=1;
        end
  end
    L(index)=laguerre(2*om(index).^ga,k,c);
    psif(:,:,k+1)=coeff.*psizero.*L;%maxmax(coeff),maxmax(psizero),maxmax(L)
end


%  See Olhede and Walden, "Noise reduction in directional signals
%  using multiple Morse wavelets", IEEE Trans. Bio. Eng., v50, 51--57.
%  The equation at the top right of page 56 is equivalent to the
%  preceding expressions. Morse wavelets are defined in the frequency  
%  domain, and so not interpolated in the time domain in the same way
%  as other continuous wavelets.


%The second family is not used for anything right now.  It's experimental.
function[X]=morsewave_second_family(fact,N,K,ga,be,om,psizero,str)

dt=1;
dom=om(2)-om(1);
a0=morseafun(ga,be,'energy');
index=(1:round(N/2));
psi0=dt.*sqrt(1./fact).*a0.*psizero;

if K>3
    error('Sorry, can only compute the first 3 members of this family right now.');
end

phi=zeros(length(om),K);
    
for k=0:K-1
    ak=morsewave_ak(k,ga,be);
    for n=0:k
        cnk=morsewave_cnk(n,k,ga,be);
        phi(index,k+1)=phi(index,k+1)+ak.*cnk.*om(index).^(n.*ga-k).*psi0(index);
    end
end

%Ensure zero mean
if iseven(N)
    phi(1,:)=0;
end

psi=phi;
for k=0:K-1
    for n=0:k-1
        bnk=morsewave_bkl(k,n,ga,be);
        psi(:,k+1)=psi(:,k+1)-bnk.*phi(:,n+1);
    end
%   morsewave_atildek(k,ga,be)
    %psi(:,k+1)=morsewave_atildek(k,ga,be).*psi(:,k+1);
    %psi(:,k+1)=psi(:,k+1)./sqrt(vsum(psi(:,k+1).^2,1))*sqrt(N);
end
psi(:,1)=2*psi(:,1)./maxmax(psi(:,1));
psi(:,2)=2*psi(:,2)./maxmax(psi(:,2));
psi(:,3)=-2*psi(:,3)./minmin(psi(:,3));
X=psi;
%X(:,1)=psi(:,1);
X(:,3)=psi(:,1)+psi(:,3);
X(:,3)=2*X(:,3)./maxmax(X(:,3));


%X=phi;

function[psif]=morsewave_second_family_alternate(fact,N,K,ga,be,om,psizero,str)

r=(2*be+1)./ga;
c=r-1;
L=0*om;
index=(1:round(N/2));

if K>2
    error('Sorry, can only compute the first two members of the edge family right now.');
end
phi=zeros(length(om),1,K);

for k=0:K-1
    ak=morsewave_ak(k,ga,be);
    for n=0:k
        cnk=morsewave_cnk(n,k,ga,be);
        phi(index,:,k+1)=phi(index,k+1)+ak.*cnk.*om(index).^(n.*ga-k).*psizero(index);
    end
end
vswap(phi,nan,0);

psif=phi;

% Older code for continuing to higher-order members of this family
% for k=0:K-1
%     for n=0:k-1
%         bnk=morsewave_bkl(k,n,ga,be);
%         psif(:,k+1)=psif(:,k+1)-bnk.*phi(:,n+1);
%     end
%    % morsewave_atildek(k,ga,be)
%     %psif(:,k+1)=morsewave_atildek(k,ga,be).*psif(:,k+1);
%     psif(:,k+1)=psif(:,k+1)./sqrt(vsum(psif(:,k+1).^2,1))*sqrt(N);
% end
% if strcmpi(str,'ban')
%     r=(2*be+1)./ga;
%     A=double((pi*ga*(2.^r)*exp(-gammaln(r))).^(1/2));
%     coeff = sqrt(2./fact)*A;
%     psif=2.*psif./coeff;
% end

if strcmpi(str,'ene')
    A=morseafun(k,ga,be,str);
    coeff = sqrt(1./fact)*A;
    psif=psif.*coeff;
end


function[cnk]=morsewave_cnk(n,k,ga,be)
if k==0&&n==0
    cnk=1;
elseif k==1&&n==0||k==0&&n==1
     cnk=be;
elseif k==1&&n==1
     cnk=-ga;        
elseif k==2&&n==0||k==0&&n==2
     cnk=be.*(be-1);
elseif k==2&&n==1||k==1&&n==2
     cnk=-(ga.*(ga-1)+2.*be.*ga); 
elseif k==2&&n==2
     cnk=ga.^2;                     
end


function[ak]=morsewave_ak(k,ga,be)
akinv=0;

for n=0:k
    for p=0:k
        cnk=morsewave_cnk(n,k,ga,be);
        cpk=morsewave_cnk(p,k,ga,be);
        ratn=frac(morseafun(ga,be,'energy'),morseafun(ga,be+n*ga-k,'energy'));
        ratp=frac(morseafun(ga,be,'energy'),morseafun(ga,be+p*ga-k,'energy'));
        akinv=akinv+cnk.*cpk.*morseproj(ga,be+n*ga-k,be+p*ga-k).*ratn.*ratp;
    end
end
ak=sqrt(1./akinv);


function[bkl]=morsewave_bkl(k,l,ga,be)
bkl=0;
for n=0:k
    for p=0:l
        cnk=morsewave_cnk(n,k,ga,be);
        cpk=morsewave_cnk(p,l,ga,be);
        ratn=frac(cnk,morseafun(ga,be+n*ga-k,'energy'));
        ratp=frac(cpk,morseafun(ga,be+p*ga-l,'energy'));
        bkl=bkl+morseproj(ga,be+n*ga-k,be+p*ga-l).*ratn.*ratp;
    end
end
bkl=bkl.*morseafun(ga,be,'energy').^2.*morsewave_ak(k,ga,be).*morsewave_ak(l,ga,be);
    

function[dnk]=morsewave_dnk(n,k,ga,be)

a1=morsewave_ak(1,ga,be);
a2=morsewave_ak(2,ga,be);
if k==0&&n==0
    dnk=1;
elseif k==1&&n==0
     dnk=-a1.*morsewave_cnk(0,1,ga,be);
elseif k==1&&n==1
     dnk=1;        
elseif k==2&&n==0
     dnk=(a2.*morsewave_cnk(0,2,ga,be).*a1-1).*morsewave_cnk(0,1,ga,be);
elseif k==2&&n==1
     dnk=-a2.*morsewave_cnk(1,2,ga,be).*a1;
elseif k==2&&n==2
     dnk=1;                     
end


%This is not yet working, and also, is really slow.
function[atildek]=morsewave_atildek(k,ga,be)

akinv=0;

for n=0:k
    for p=0:k
        dnk=morsewave_dnk(n,k,ga,be);
        dpk=morsewave_dnk(p,k,ga,be);
        %ratn=frac(morseafun(ga,be,'energy'),morseafun(ga,be+n*ga-k,'energy'));
        %ratp=frac(morseafun(ga,be,'energy'),morseafun(ga,be+p*ga-k,'energy'));
        akinv=akinv+dnk.*dpk.*morsewave_bkl(n,p,ga,be);%.*ratn.*ratp;
    end
end
atildek=sqrt(1./akinv);

function[y]=laguerre(x,k,c)
%LAGUERRE Generalized Laguerre polynomials
%
%   Y=LAGUERRE(X,K,C) where X is a column vector returns the
%   generalized Laguerre polynomials specified by parameters K and C.
%  
%   LAGUERRE is used in the computation of the generalized Morse
%   wavelets and uses the expression given by Olhede and Walden (2002),
%  "Generalized Morse Wavelets", Section III D. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 F. Rekibi and J. M. Lilly 
%                         --- type 'help jlab_license' for details  


y=0*x;
for m=0:k
   %Log of gamma function much better ... trick from Maltab's ``beta'' 
   %fact=gamma(k+c+1)./(gamma(c+m+1).*gamma(k-m+1));  
   fact=exp(gammaln(k+c+1)-gammaln(c+m+1)-gammaln(k-m+1));
   y=y+((-1).^m).*fact.*(x.^m)./gamma(m+1);
end



function[]=morsewave_test
morsewave_test_unitenergy
morsewave_test_centering
morsewave_test_scorer
morsewave_test_cauchy
morsewave_test_gaussian
morsewave_test_jdawson
morsewave_test_admiss
morsewave_test_fft



function[]=morsewave_test_unitenergy
fs=2*pi./logspace(log10(5),log10(40))'; 
N=1023;
w=morsewave(N,2,2,4,fs,'energy');
energy=vsum(abs(w(:,:,1)).^2,1);
reporttest('MORSEWAVE unit energy for unit sample rate, K=1',aresame(energy,1+0*energy,1e-4))
energy=vsum(abs(w(:,:,2)).^2,1);
reporttest('MORSEWAVE unit energy for unit sample rate, K=2',aresame(energy,1+0*energy,1e-4))

% w=morsewave(N,2,4,fs,'energy','edge');
% energy=vsum(abs(w(:,:,2)).^2,1);
% reporttest('MORSEWAVE unit energy for edge wavelet',maxmax(energy-1)<1e-4)
% energy=vsum(w(:,:,2).*conj(w(:,:,1)),1);
% reporttest('MORSEWAVE orthogonality for edge wavelet',maxmax(abs(energy))<1e-6)


function[]=morsewave_test_centering
fs=2*pi./logspace(log10(5),log10(40))'; 
N=1023;
w=morsewave(N,2,4,fs);
bool=0*fs;
for i=1:size(w,2)
   bool(i)=max(abs(w(:,i)))==abs(w(N/2+1/2,i));
end
reporttest('MORSEWAVE centered for odd N',all(bool))

N=1024;
w=morsewave(N,2,4,fs);
bool=0*fs;
for i=1:size(w,2)
   bool(i)=max(abs(w(:,i)))==abs(w(N/2,i)) || max(abs(w(:,i)))==abs(w(N/2+1,i));
end
reporttest('MORSEWAVE centered for even N',all(bool))


function[]=morsewave_test_scorer

dt=1;
t=(-50:dt:50)';
psi1=morsewave(length(t),1,3,0,morsefreq(3,1)./dt,'bandpass');
c=3.^(1/3);
psi2=(1./c).*scorer(sqrt(-1)*t./c,1000);
%figure,uvplot(psi1),
%figure,uvplot(psi2)

err=vsum(abs(psi1-psi2).^2,1)./vsum(abs(psi1).^2,1);
reporttest('MORSEWAVE for GAMMA=3 wavelet versus Scorer function expression',err<1e-1)

function[hi]=scorer(z,n)
%This is not a very good way to compute the Scorer functions.
%Need to set n really high to have to integal behave nicely.
%Only used for testing purposes.

z=z(:);
u=linspace(0,10,n);
du=u(2)-u(1);
umat=vrep(u,size(z,1),1);

%aiz=airy(0,z);
%biz=airy(2,z);
%gi=frac(1,pi).*vsum(sin(frac(tmat.^3,3)+oprod(z,t)),2).*dt;
hi=frac(1,pi).*vsum(exp(-frac(umat.^3,3)+oprod(z,u')),2).*du;


function[]=morsewave_test_cauchy


dt=.01;
t=(-25:dt:25)';
psi1=morsewave(length(t),1,1,0,morsefreq(1,1).*dt,'bandpass')./dt;
psi2=frac(1,pi).*frac(1,1-sqrt(-1)*t);

err=vsum(abs(psi1-psi2).^2,1)./vsum(abs(psi1).^2,1);
reporttest('MORSEWAVE for GAMMA=1 wavelet versus Cauchy function expression',err<1e-1)


function[]=morsewave_test_gaussian

dt=.1;
t=(-50:dt:50)';
psi1=morsewave(length(t),1,2,0,morsefreq(2,1).*dt,'bandpass')./dt;
psi2=frac(1,2*sqrt(pi)).*(exp(-(t/2).^2)+sqrt(-1)*jdawson(t/2)*frac(2,sqrt(pi)));
err=vsum(abs(psi1-psi2).^2,1)./vsum(abs(psi1).^2,1);
reporttest('MORSEWAVE for GAMMA=2 wavelet versus Gaussian function expression',err<1e-1)

function[]=morsewave_test_jdawson
dt=0.01;
t=(-15:dt:15)';

n=5;
herm=hermpoly(t(:)/2,n+1);
herm=herm(:,2:end);
g=exp(-frac(t.^2,4));

[psi1,psi2]=vzeros(length(t),5);
for k=1:5
    dk=jdawson(t/2,k);
    coeffk=frac(1,4*sqrt(pi)).*morseafun(2,k).*frac(sqrt(-1),2).^k;
    tic
    psi1(:,k)=morsewave(length(t),1,2,k,morsefreq(2,k).*dt,'bandpass')./dt;
    toc
    tic
    psi2(:,k)=coeffk*(g.*herm(:,k)+sqrt(-1)*(-1).^k.*dk*frac(2,sqrt(pi)));
    toc
    err=vsum(abs(psi1(:,k)-psi2(:,k)).^2,1)./vsum(abs(psi1(:,k)).^2,1);
    reporttest(['MORSEWAVE for GAMMA=2 derivatives matches Dawson expression for n=' int2str(k)],err<1e-3)
end


function[]=morsewave_test_admiss
ga1=(1:1:11);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);

N=1000;
psi2om=zeros(N,length(ga));
dt=1/10;
om=(0:N-1)'./N;om(1)=1e-4;
for i=1:length(ga)
    [psii,psifi]=morsewave(N,1,ga(i),be(i),morsefreq(ga(i),be(i)).*dt,'bandpass');
    psi2om(:,i)=(psifi).^2./om;
end
cpsi1=vsum(psi2om,1).*(1/N);
cpsi1=reshape(cpsi1,length(be1),length(ga1));
cpsi2=morseafun(ga,be).^2.*frac(1,ga.*2.^(2*be./ga)).*gamma(frac(2*be,ga));
cpsi2=reshape(cpsi2,length(be1),length(ga1));
        
%ga1=(1:1:11);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);

N=1000;
psi2om=zeros(N,length(ga));
dt=1/10;
om=(0:N-1)'./N;om(1)=1e-4;
for i=1:length(ga)
    [psii,psifi]=morsewave(N,1,ga(i),be(i),morsefreq(ga(i),be(i)).*dt,'bandpass');
    psi2om(:,i)=(psifi).^2./om;
end
cpsi1=vsum(psi2om,1).*(1/N);
cpsi1=reshape(cpsi1,length(be1),length(ga1));
cpsi2=morseafun(ga,be).^2.*frac(1,ga.*2.^(2*be./ga)).*gamma(frac(2*be,ga));
cpsi2=reshape(cpsi2,length(be1),length(ga1));
ga1=(1:1:11);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);

N=1000;
psi2om=zeros(N,length(ga));
dt=1/10;
om=(0:N-1)'./N;om(1)=1e-4;
for i=1:length(ga)
    [psii,psifi]=morsewave(N,1,ga(i),be(i),morsefreq(ga(i),be(i)).*dt,'bandpass');
    psi2om(:,i)=(psifi).^2./om;
end
cpsi1=vsum(psi2om,1).*(1/N);
cpsi1=reshape(cpsi1,length(be1),length(ga1));
cpsi2=morseafun(ga,be).^2.*frac(1,ga.*2.^(2*be./ga)).*gamma(frac(2*be,ga));
cpsi2=reshape(cpsi2,length(be1),length(ga1));
[p,dt,dom]=morsebox(ga,be);

reporttest('MORSEWAVE admissibility matches analytic expression',aresame(cpsi1,cpsi2,1e-2));




function[]=morsewave_test_fft
fs=1./logspace(log10(5),log10(40))'; 
N=1023;
[psi,Psi]=morsewave(N,2,2,4,fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N odd, FS positive, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morsewave(N+1,2,2,4,fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N even, FS positive, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morsewave(N,2,1,4,-fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N odd, FS negative, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morsewave(N+1,2,2,4,-fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N even, FS negative, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))

% Figure for second family, not currently working 
% %function[]=morsewave_figure2
% fs=2*pi*.05;N=1000;
% 
% 
% [x2,X2]=morsewave(N,3,2,4,fs,'energy','first');
% [y2,Y2]=morsewave(N,3,2,4,fs,'energy','edge');
% 
% [x3,X3]=morsewave(N,3,3,4,fs,'energy','first');
% [y3,Y3]=morsewave(N,3,3,4,fs,'energy','edge');
% 
% [x4,X4]=morsewave(N,3,4,4,fs,'energy','first');
% [y4,Y4]=morsewave(N,3,4,4,fs,'energy','edge');
% 
% figure
% subplot(231),
% plot((1:length(X2))/N,X2./maxmax(X2(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
% title('First Family | \gamma=2, \beta=4'),ylim([-1.1 1.1])
% subplot(232),
% plot((1:length(X2))/N,X3./maxmax(X3(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
% title('First Family | \gamma=3, \beta=4'),ylim([-1.1 1.1])
% subplot(233),
% plot((1:length(X2))/N,X4./maxmax(X4(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
% title('First Family | \gamma=4, \beta=4'),ylim([-1.1 1.1])
% subplot(234),
% plot((1:length(X2))/N,Y2./maxmax(Y2(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
% title('Second Family | \gamma=2, \beta=4'),ylim([-1.1 1.1])
% subplot(235),
% plot((1:length(X2))/N,Y3./maxmax(Y3(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
% title('Second Family | \gamma=3, \beta=4'),ylim([-1.1 1.1])
% subplot(236),
% plot((1:length(X2))/N,Y4./maxmax(Y4(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
% title('Second Family | \gamma=4, \beta=4'),ylim([-1.1 1.1])
% 
% packfig(2,3,'columns');
