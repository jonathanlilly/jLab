function[varargout]=transmaxdist(varargin)
%TRANSMAXDIST  Distributions of wavelet transform maxima in noise.
%
%   This function is part of 'element analysis' described in Lilly (2017), 
%   "Element analysis: a wavelet-based method for analyzing time-localized
%   events in noisy time series", submitted.  Available at www.jmlilly.net.
%  
%   [COUNT,BINS]=TRANSMAXDIST(GAMMA,BETA,ALPHA,FS,R,N,M) returns the 
%   histogram of wavelet transform maxima magnitudes, for a length N time 
%   series having spectral slope -2*ALPHA transformed at frequencies FS 
%   using a (GAMMA,BETA) wavelet, based on a simulation having N*M points.
%
%   Here GAMMA, BETA, ALPHA, and R are all scalars, or are all arrays of 
%   the same length as FS.  FS is a frequency array computed by MORSESPACE.
%
%   R is the ratio between each frequency FS and the next.  This will be 
%   constant and greater than one when FS is computed by MORSESPACE.  As 
%   as described in Appendix C of Lilly (2017), R=FS(n)./FS(n+1) for all n. 
%
%   COUNT is the number of transform maxima observed at each frequency in
%   the magnitude bins BINS.  COUNT is a LENGTH(BINS) x LENGTH(FS) matrix.
%
%   Transform maxima values, as output in BINS, are normalized such that
%   the expected squared magnitude of the wavelet transform of noise occurs
%   at unity.  BINS thus corresponds to the normalized event magnitude.  
%
%   TRANSMAXDIST works by simulating a vector whose statistical properties 
%   mimic those of the wavelet transform and the four adjacted points, thus 
%   avoiding the need to explicitly compute the transform.  The choice of
%   e.g. M=1000 simulates a transform 1000 times as long as time series of 
%   interest, which itself is of length N. 
%
%   [COUNT,BINS,RATE]=TRANSMAXDIST(...) also returns the RATE, the
%   normalized reversed cumulative density function.  RATE gives the 
%   expected number of transform maxima occuring in a time series of length
%   N having a magnitude greater than the corresponding bin value.
%
%   [COUNT,BINS,RATE,SIGMA]=TRANSMAXDIST(...) also returns the theoretical 
%   covariance matrix SIGMA from which the Monte Carlo simulations are 
%   constructed.  SIGMA is an array of length 5 x 5 x LENGTH(FS). 
%
%   Note that if the covariance matrix is not positive definite, as can 
%   happen due to numerical complications for extreme BETA and GAMMA 
%   choices, then COUNT and RATE will both consist entirely of NaNs.
%   _______________________________________________________________________
%
%   Additional options
%
%   TRANSMAXDIST(...,BINS) alternately uses BINS for the bin centers 
%   instead of the default choice, which is set to LINSPACE(0,6,400)'.
%
%   By default, TRANSMAXDIST performs a simulation for each of the scale
%   frequencies in FS.  TRANSMAXDIST(...,'extrapolate') instead computes
%   the distribution only for the highest scale frequency, then 
%   extrapolates these values to all other scale frequencies with a scaling
%   law.  GAMMA, BETA, ALPHA, and R must all be scalars in this case.
%  
%   For details, see Lilly (2017).
%
%   See also MAXPROPS, TRANSMAX, ISOMAX, MAX2EDDY.
%
%   'transmaxdist --t' runs some tests.
%
%   Usage: [count,bins]=transmaxdist(ga,be,al,fs,r,N,M); 
%          [count,bins,rate,sigma]=transmaxdist(ga,be,al,fs,r,N,M);
%          [count,bins,rate,sigma]=transmaxdist(ga,be,al,fs,r,N,M,'extrap');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    transmaxdist_test,return
end

%parstr='serial';
%
%   TRANSMAXDIST(...,'parallel') uses a PARFOR loop in the computation of 
%   the covariance matrix SIGMA.  The default behavior is 'serial'.

str='all';
normstr='band';
for i=1:3
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'all')||strcmpi(varargin{end}(1:3),'ext')
            str=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'ser')||strcmpi(varargin{end}(1:3),'par')
            parstr=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'ban')||strcmpi(varargin{end}(1:3),'ene')
            normstr=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

gamma=varargin{1};
beta=varargin{2};
alpha=varargin{3};
fs=varargin{4};
r=varargin{5};
N=varargin{6};
M=varargin{7};
if length(varargin)==7
    bins=linspace(0,6,200)';
else
    bins=varargin{8}(:);
end

s=morsefreq(gamma,beta)./fs;

%Some error checking
if ~isscalar(gamma)
    lg=length(gamma);
    lb=length(beta);
    la=length(alpha);
    ls=length(s);
    lr=length(r);
    if ~(lg==lb&&lg==la&&lg==ls&&lg==lr)
        error('TRANSMAXDIST was expecting the first five input arguments to all be the same length.')
    end
    if strcmpi(str(1:3),'ext')&&(length(gamma)>1)
        error('Sorry, ALPHA, BETA, GAMMA, FS, and R must be scalars for ''extrapolate'' option.')
    end
else
    if length(s)>1
        ro=s(2:end)./s(1:end-1);
        if ~allall(abs(ro-r)<1e8)
            disp('Input value of R does not match that computed from FS.');
        end
    end
end

arrayify(s,r,gamma,beta,alpha);
Sigma=zeros(5,5,length(s));

%s,r,gamma,beta,alpha
%if strcmpi(parstr(1:3),'ser')
    for i=1:length(s)
        Sigma(:,:,i)=sigmamat1(s(i),r(i),gamma(i),beta(i),alpha(i),normstr);
    end
%elseif strcmpi(parstr(1:3),'par')
%    parfor i=1:length(s)
%        Sigma(:,:,i)=sigmamat1(s(i),r(i),gamma(i),beta(i),alpha(i),normstr);
%    end
%end

xeps=(randn(5,N*M)+1i*randn(5,N*M))./sqrt(2);%std(xeps(:))
n=zeros(length(bins),length(s));

if strcmpi(str(1:3),'all')
    for k=1:length(s)
        disp(['TRANSMAXDIST performing simulation ' int2str(k) ' of ' int2str(length(s)) '.'])
        try
           L=chol(Sigma(:,:,k),'lower');
        catch
           L=[];
        end
        
        if ~isempty(L)
            wsim=L*xeps;
            
            bool=(abs(wsim(1,:))>abs(wsim(2,:)));
            bool=bool&(abs(wsim(1,:))>abs(wsim(3,:)));
            bool=bool&(abs(wsim(1,:))>abs(wsim(4,:)));
            bool=bool&(abs(wsim(1,:))>abs(wsim(5,:)));
            
            wmax=wsim(1,bool);
            nk=hist(abs(wmax),bins);
            n(:,k)=nk';
        else
            n(:,k)=nan*ones(size(n(:,k)));
        end
    end
elseif strcmpi(str(1:3),'ext')
    L=chol(Sigma(:,:,1),'lower');
    wsim=L*xeps;
    
    bool=(abs(wsim(1,:))>abs(wsim(2,:)));
    bool=bool&(abs(wsim(1,:))>abs(wsim(3,:)));
    bool=bool&(abs(wsim(1,:))>abs(wsim(4,:)));
    bool=bool&(abs(wsim(1,:))>abs(wsim(5,:)));
    
    wmax=wsim(1,bool);
    nk=hist(abs(wmax),bins);
    n(:,1)=nk';
    for k=2:length(s)
        n(:,k)=n(:,1)*frac(s(1),s(k));
    end
end

n=n./M;
rate=cumsum(n,1,'reverse');

%rate=n;
%for i=1:size(n,2)
    %L=2*sqrt(2)*sqrt(gamma(i).*beta(i))./fs(i);
 %   rate(:,i)=cumsum(n(:,i),1,'reverse');
%end

%     for i=1:5
%         for j=1:5
%             Sigmahat(i,j)=vmean(wsim(i,:).*conj(wsim(j,:)),2);
%         end
%     end
%     Sigmahat=Sigmahat./vmean(squared(wsim(1,:)),2);

varargout{1}=n;
varargout{2}=bins;
varargout{3}=rate;
varargout{4}=Sigma;
%varargout{4}=Sigmahat;

function[Sigma]=sigmamat1(s,r,gamma,beta,alpha,normstr)

p={gamma,beta,alpha,normstr};
Sigma=zeros(5,5);

Sigma(1,:)=[xi(0,s,1,p) xi(1,s,1,p)  xi(-1,s,1,p)  xi(0,s,r,p)   xi(0,s,1/r,p)     ];
Sigma(2,:)=[nan         xi(0,s,1,p)  xi(-2,s,1,p)  xi(-1,s,r,p)  xi(-1,s,1/r,p)    ];
Sigma(3,:)=[nan         nan          xi(0,s,1,p)   xi(1,s,r,p)   xi(1,s,1/r,p)     ];
Sigma(4,:)=[nan         nan          nan           xi(0,r*s,1,p) xi(0,r*s,1/r.^2,p)];
Sigma(5,:)=[nan         nan          nan           nan           xi(0,s/r,1,p)     ];

for j=1:5
    for i=(j+1):5
          Sigma(i,j)=conj(Sigma(j,i));
    end
end

[m0,ffun]=morsemom(-2*alpha,gamma,beta);
if strcmpi(normstr(1:3),'ene')
    Sigma=frac(Sigma,ffun.*s.^(2*alpha));
elseif strcmpi(normstr(1:3),'ban')
    Sigma=frac(Sigma,ffun*s.^(2*alpha-1));
end
    
function[x]=xi(tau,s,r,p)

gamma=p{1};
beta=p{2};
alpha=p{3};
normstr=p{4};

rtilde=(1+r.^gamma).^(1./gamma);

%Implicitly set A^2=1 as I will shortly divide by it
%fact1=frac(morseafun(gamma,beta).^2,morseafun(gamma,2*beta-2*alpha));
%fact2=frac((r.^beta).*(s.^(2*alpha-1)),rtilde.^(2*beta-2*alpha+1));

%Note:  use the input normalization in the numerator, and the *amplitude*
%normalization in the denominator.  All the latter does is cancel the 
%coefficient coming from MORSEXPAND, which assumes the amplitude normalization
fact1=frac(morseafun(gamma,beta,normstr).^2,morseafun(gamma,2*beta-2*alpha));
if strcmpi(normstr(1:3),'ene')
    fact2=frac((r.^(beta+1/2)).*(s.^(2*alpha)),rtilde.^(2*beta-2*alpha+1));
elseif strcmpi(normstr(1:3),'ban')
    fact2=frac((r.^beta).*(s.^(2*alpha-1)),rtilde.^(2*beta-2*alpha+1));
end
    
x=fact1.*fact2.*conj(morsexpand(tau./(s.*rtilde),gamma,2*beta-2*alpha,morsefreq(gamma,2*beta-2*alpha)));
%x=fact1.*fact2.*conj(morsexpand(tau./(s.*rtilde),gamma,2*beta-2*alpha,morsefreq(gamma,2*beta-2*alpha),'cumulant'));

%fact1,fact2
%Equivalent
%fact2=frac((r.^beta).*(s.^(2*alpha)),rtilde.^(2*beta-2*alpha));
%fs=morsefreq(gamma,2*beta-2*alpha)./(s.*rtilde);
%x=fact1.*fact2.*conj(morsexpand(tau,gamma,2*beta-2*alpha,fs));


function[]=transmaxdist_test

for k=1:2
    tic
    ga=2;be=2;
    switch k
        case 1
            alpha=0;
        case 2
            alpha=1;
    end
    %fs=morsespace(ga,be,100);
    fs=morsespace(ga,be,{0.01,pi},2*pi/100);
    N=1e6;
    rng(1);
    x=randn(N,1);
    if alpha ==1
        x=cumsum(x);
        x=x./std(x);
    end
    w=wavetrans(x,{ga,be,fs(1:3),'bandpass'});
    
    clear xvec
    xvec(:,1)=w(:,2);
    xvec(:,2)=w([2:end 1],2);
    xvec(:,3)=w([end 1:end-1],2);
    xvec(:,4)=w(:,3);
    xvec(:,5)=w(:,1);
    
    for i=1:5
        for j=1:5
            Sigmahat(i,j)=vmean(xvec(:,i).*conj(xvec(:,j)),1);
        end
    end
    Sigmahat=Sigmahat./vmean(squared(xvec(:,1)),1);
    
    [count,bins,rate,Sigma]=transmaxdist(ga,be,alpha,fs(2),fs(1)./fs(2),N,1);
    switch k
        case 1
            bool=maxmax(abs(Sigma-Sigmahat)./abs(Sigma))<1/100;
            reporttest('TRANSMAXDIST simulated and theoretical covariance matrices agree to within 1%, white noise case',bool)
        case 2
            bool=maxmax(abs(Sigma-Sigmahat)./abs(Sigma))<4/100;
            reporttest('TRANSMAXDIST simulated and theoretical covariance matrices agree to within 4%, red noise case',bool)
    end
    
    %[index,ww]=transmax(fs(1:3),w);
    %n=hist(abs(ww)./vstd(w(:,2),1),bins)';
    %figure,plot(bins,[count n]) %That looks great!
end

function[]=transmaxdist_other


%[n,sigma]=transmaxdist(ga,be,0,[fs(2) fs(3)],N/10);
N=1e8;
tic;[n,x,Sigma]=transmaxdist(ga,be,0,fs,N);toc
for i=1:size(n,2)
    ntilde(:,i)=n(:,i)./frac(N,2*sqrt(2)*sqrt(ga*be)./fs(i));
    %ntilde(:,i)=n(:,i)./sum(n(:,i));
end


tic
N=1e8;
x=randn(N,1);
w=wavetrans(x,{ga,be,fs(1:3)});
toc

tic; [ii,jj,ww,ff]=transmax(fs(1:3),w);toc
bins=linspace(0,6,200)';
wwtilde=ww./sqrt(vmean(squared(w(:,2)),1));
nn=hist(abs(wwtilde),bins);
nn=nn./frac(N,2*sqrt(2)*sqrt(ga*be)./fs(2));



