function[varargout]=mspec(varargin)
% MSPEC  Multitaper power and cross spectra.
%   _______________________________________________________________________
%
%   *|* mspec.png --- Figure illustrating the multitaper spectrum. 
%   Type 'jhelp mspec' to view this image. *|*
%   _______________________________________________________________________
%
%   MSPEC implements spectral and cross-spectral analysis using the 
%   multitaper method for real or complex-valued data. 
%
%   MSPEC is to be run after calling SLEPTAP to compute the multitapers.
%   _______________________________________________________________________
%
%   One-sided power spectrum
%
%   [F,S]=MSPEC(X,PSI) returns the one-sided power spectrum of X at 
%   positive frequencies using data tapers PSI. 
%  
%   Input:   
%       X  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
% 
%   Output:
%       F  --  [M/2] nonnegative *radian* frequencies
%       S  --  [M/2] x N one-sided power spectrum matrix
%                
%   In the above, [M/2] means: M/2 if M is even, and (M+1)/2 if M is odd. 
%   See FOURIER for the calculation of the Fourier frequencies.
%
%   The spectra matrices are averages over the K "eigenspectra" computed 
%   with each of the K tapers, as discussed in Park et al. JGR 1987.
%
%   The one-sided spectrum S is normalized such that its sum over all 
%   frequencies F, (1/2/pi)*SUM(S,1)*DF where DF is the frequency 
%   increment, approximates the signal variance; see discussion below.
%
%   MSPEC(...,'detrend') detrends the data before computing the spectra.
%   ______________________________________________________________________
%  
%   Cross-spectra of real-valued data
%   
%   MSPEC can be used to compute the cross-spectrum of two real-valued
%   time series or sets of time series.
%
%   [F,SXX,SYY,SXY]=MSPEC(X,Y,PSI);   --- For cross-spectra
%        
%   Input:
%       X  --  M x N matrix containing N length M time series
%       Y  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
% 
%   Output:
%       F  --  M/2 nonnegative frequencies
%     SXX  --  [M/2] x N one-sided spectra of X
%     SYY  --  [M/2] x N one-sided spectra of Y
%     SXY  --  [M/2] x N one-sided cross spectra of X and Y 
%
%   See TWOSPECPLOT for plotting SXX and SYY simultaneously.
%   ______________________________________________________________________
%  
%   Rotary spectra of complex-valued data
%
%   MSPEC can also the so called "rotary spectra" of complex-valued 
%   time series or sets of time series.
%
%   [F,SPP,SNN,SPN]=MSPEC(Z,PSI);   --- For rotary spectra of Z=X+iY
%
%   Input:   
%       Z  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
%   
%   Output:   
%        F  --  M/2 nonnegative frequencies
%      SPP  --  [M/2] x N positively rotating power spectrum matrix
%      SNN  --  [M/2] x N negatively rotating power spectrum matrix  
%      SPN  --  [M/2] x N rotary cross spectral matrix      
%
%   Note that the rotary spectra are defined such that SXX+SYY=SPP+SNN.
%  
%   The rotary spectra SPP and SNN are normalized such that the sum of SPP
%   over all frequencies plus that of SNN approximates the variance of Z. 
%
%   See TWOSPECPLOT for plotting SPP and SNN simultaneously.
%   ______________________________________________________________________
%  
%   Sample rate
%
%   [F,S]=MSPEC(DT,...) specifies the sample interval to be used in the
%   calculation of the frequency array F. DT defaults to unity.
%
%   Spectral values depend linearly upon the sample rate in order that the 
%   integral of the spectra over frequency approximate the variance.
%   ______________________________________________________________________
%
%   Periodogram
%
%   MPSEC can be used to form the naive spectral estimator, known as the
%   periodogram. Although this is not generally a good way to estimate the
%   spectrum, it can be useful as a comparision.
%
%   MSPEC(X,[]) or MSPEC(X,Y,[]) with PSI empty uses the default, or boxcar
%   taper, normalized to unit energy. This returns the periodogram.  
%   ______________________________________________________________________
%  
%   Normalizations
%
%   By default, MSPEC uses *radian* frequency as in cos(f t).  Optionally
%   MSPEC(,...,'cyclic') will use *cyclic* frequency, as in cos(2 pi f t).
%
%   MSPEC is normalized to approximately recover the time series variance. 
%   For the MSPEC periodogram, this recovery is exact, although the 
%   expressions are complicated somewhat by the use of one-sided spectra.  
%
%   Real-valued data
%
%   [F,S]=MSPEC(DT,X,[]) where X is a real-valued time series of length M
%   recovers the variance of X, STD(X,1).^2, as follows:
%
%     (1/2/pi)*(F(2)-F(1))*SUM(S(2:end))               -- M odd
%     (1/2/pi)*(F(2)-F(1))*(SUM(S(2:end-1))+S(end)/2)  -- M even
%
%   Note that the zero frequency is omitted in the summation, and for even 
%   time series length, the power at the Nyquist S(end) must be divided by 
%   two to avoid double-counting by the one-sided spectrum.  The "1" in the 
%   argument of STD forces STD to use an N rather than N-1 normalization. 
%
%   Complex-valued data
%
%   [F,SPP,SNN]=MSPEC(DT,Z,[]) where Z is a complex-valued time series of 
%   length M recovers the variance of Z, STD(Z,1).^2, as follows:
%
%     (1/2/pi)*(F(2)-F(1))*(SUM(SPP(2:end))+SUM(SNN(2:end)))   -- M odd
%     (1/2/pi)*(F(2)-F(1))*(SUM(SPP(2:end))+SUM(SNN(2:end-1))) -- M even
%
%   Again the modification for even M prevents the power at the Nyquist
%   from being double-counted.  This modification is necessary because the 
%   negative rotary spectrum duplicates the Nyquist when M is even.  
%   ______________________________________________________________________
%   
%   Cross-spectra of complex-valued data
%
%   To compute the cross-spectra of two complex-valued time series or sets 
%   of time series Z1 and Z2, run MSPEC repeatedly.
%
%   [F,SP1P1,SP2P2,SP1P2]=MSPEC(Z1,Z2,PSI);  
%   [F,SN1N1,SN2N2,SN1N2]=MSPEC(CONJ(Z1),CONJ(Z2),PSI);  
%
%   The first call returns the spectra and cross-spectra of Z1 and Z2 at
%   positive frequencies, while the second returns their spectra and cross-
%   spectra at negative frequencies.  Finally
%
%   [F,SP1P1,SN2N2,SP1N2]=MSPEC(Z1,CONJ(Z2),PSI);  
%
%   returns the cross-spectra between the positive rotary components of Z1 
%   and the negative components of Z2.
%   ______________________________________________________________________
%
%   Adaptive spectra
%
%   MSPEC(...,LAMBDA,'adaptive'), where LAMBDA contains the eigenvalues of
%   the tapers as computed by SLEPTAP, alternately uses the "adaptive"
%   multitaper method of Thomson (1982).
% 
%   This implementation follows that of Park et al. (1987a), JGR.
%
%   For cross-spectra or for rotary spectra, the weights appearing in the
%   adaptive spectra are derived for the total spectrum of each signal 
%   compoment, that is for SXX+SYY or SPP+SNN as appropriate.  Then the
%   separate spectra and co-spectra are computed using identical weights.
%   ______________________________________________________________________
%  
%   Cell array input / output
%
%   MSPEC generates cell array output given cell array input.
%
%   Let's say one has P different time series, X1, X2,..., XP.  Put these 
%   into a cell array X{1}=X1, X{2}=X2, ..., X{P}=XP, and then use
%   '[psi,lambda]=sleptap(cellength(x))' to make a cell array of tapers.
%
%   [F,S]=MSPEC(X,PSI) then returns cell arrays F and S corresponding 
%   to the Fourier frequencies and spectra of the P arrays.  
%
%   The other argument forms given above also work.  In particular, 
%   specifiying the sample time through MPSEC(DT,...) works, with DT either
%   a scalar or an array of the same length as the cell array X.
%
%   The spectra can then be plotted with CELLPLOT(F,S), or TWOSPECPLOT for
%   a pair of output spectra.
%   ______________________________________________________________________
%
%   Parallelization
%
%   MSPEC(..., 'parallel') when the input fields X, X and Y, or Z are cell
%   arrays, parallelizes the spectral estimation by looping over the cells
%   with a PARFOR loop.  This requires Matlab's Parallel Computing Toolbox.
%   ______________________________________________________________________
% 
%   Example 
%
%   The example at the top of this help file shows clockwise (left) and 
%   counterclockwise (right) rotary spectra from moored current meter  
%   measurements of the ocean currents in the Labrador Sea.
%   
%   The periodogram is in gray, and blue and red are multitaper spectra 
%   with P=4 and P=32, respectively.  The local Coriolis frequency is 
%   marked with a dashed line.  Tidal and inertial peaks are apparent.
%
%   The main point of this figure is to show that increasing P increases
%   the degree of frequency-domain smoothing.   
%   ______________________________________________________________________
%
%   'mspec --t' runs some tests.
%   'mspec --f' generates the above sample figure from Bravo mooring data.
%
%   See also: SLEPTAP, HERMFUN, MTRANS, MSVD, TWOSPECPLOT.
%
%   Usage   [f,s]=mspec(x,psi);    
%           [f,s]=mspec(dt,x,psi);     
%           [f,spp,snn,spn]=mspec(z,psi);     
%           [f,sxx,syy,sxy]=mspec(x,y,psi);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details        
  
%   The Cartesian and rotary spectra are then related by a unitary
%   transformation.  For details see
%
%       Lilly and Olhede (2010).  Bivariate instantaneous frequency and
%           bandwidth.  IEEE Trans. Sig. Proc.

if strcmpi(varargin{1},'--t')
    mspec_test;return
end
if strcmpi(varargin{1},'--f')
    type makefigs_mspec;
    makefigs_mspec;
    return
end

deltat=1;
if isscalar(varargin{1})
  deltat=varargin{1};
  varargin=varargin(2:end);
else
  if nargin>2
      bool1=~iscell(varargin{2})&&(numel(varargin{1})==size(varargin{2},2));
      bool2= iscell(varargin{2})&&(numel(varargin{1})==length(varargin{2}));
      if bool1||bool2
          deltat=varargin{1};
          deltat=deltat(:)';
          varargin=varargin(2:end);
      end
  end
end
  

%Sort out input arguments
lambda=1;  %This means use the average multitaper spectrum
detrendstr='nodetrend';
normstr='rad';
cores='serial';

for i=1:4
    if ischar(varargin{end})
        if ~isempty(strfind(varargin{end},'ada'))
            lambda=varargin{end-1};
            varargin=varargin(1:end-2);
            if length(lambda)~=numel(lambda)
                error('Looks like you forgot to input LAMBDA with the adaptive algorithm.')
            end
        elseif ~isempty(strfind(varargin{end},'det'))||~isempty(strfind(varargin{end},'nod'))
             detrendstr=varargin{end};
             varargin=varargin(1:end-1);
        elseif strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
             cores=varargin{end};
             varargin=varargin(1:end-1);
        else
             normstr=varargin{end};
             varargin=varargin(1:end-1);
        end
    end
end

if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
end


x=varargin{1};
if length(varargin)==1
    error('Taper not specified.')
end
psi=varargin{end};
na=length(varargin);

y=[];
if na==3
    y=varargin{2};
end

if iscell(x)||iscell(y)
    %All of this is just to handle cell array input
    if isempty(y)
        y=cell(size(x));
    elseif isempty(x)
        x=cell(size(y));
    end
    if ~iscell(psi);
        psio=psi;
        psi=cell(size(x));
        for i=1:length(x)
            psi{i}=psio;
        end
    end
    if size(lambda,2)==1
        lambda=vrep(lambda,length(x),2);
    end
    if length(deltat)==1
        deltat=deltat+zeros(size(x));
    end
    
    if strcmpi(cores(1:3),'ser')
        for i=1:length(x)
            disp(['SLEPTAP computing spectra for time series #' int2str(i) ' of ' int2str(length(x)) '.'])
            %iscell(deltat),iscell(x),iscell(y),iscell(psi)
            %vsize(deltat(i),x{i},y{i},psi{i},lambda(:,i))
            cellout{i}=mspec_one(deltat(i),x{i},y{i},psi{i},lambda(:,i),detrendstr,normstr);
        end
    elseif strcmpi(cores(1:3),'par')
        %Exactly the same but with a parfor
        parfor i=1:length(x)
            disp(['SLEPTAP computing spectra for time series #' int2str(i) ' of ' int2str(length(x)) '.'])
            %iscell(deltat),iscell(x),iscell(y),iscell(psi)
            %vsize(deltat(i),x{i},y{i},psi{i},lambda(:,i))
            cellout{i}=mspec_one(deltat(i),x{i},y{i},psi{i},lambda(:,i),detrendstr,normstr);
        end
    end
    
    for i=1:length(x)
        for j=1:length(cellout{i})
            varargout{j}{i,1}=cellout{i}{j};
        end
    end
else
    varargout=mspec_one(deltat,x,y,psi,lambda,detrendstr,normstr);
end


function[cellout]=mspec_one(dt,x,y,psi,lambda,detrendstr,normstr)

if ~isscalar(dt)
    dt=vrep(dt,length(fourier(size(x,1))),1);
end

if ~isreal(x)&&isempty(y)
    y=conj(x);
end

if strcmpi(detrendstr(1:3),'det')
    x=detrend(x);
    if ~isempty(y)
        y=detrend(y);
    end
end

if isempty(psi)
    psi=frac(1,sqrt(size(x,1)))+zeros(size(x(:,1)));
end

%In real cases, multiply by sqrt(2) to get one-sided spectrum
if isempty(y)        %One real-valued
     [f,mmatx]=mtrans(sqrt(2)*x,psi);
else
     if ~isreal(x)   %Two complex-valued
        [f,mmatx,mmaty]=mtrans(x,y,psi);
     else            %Two real-valued
        [f,mmatx,mmaty]=mtrans(sqrt(2)*x,sqrt(2)*y,psi); 
     end
end

if isempty(y) %One time series
     if lambda==1
         cellout{2}=avgspec(mmatx,mmatx).*dt;
     else
         var=squared(vstd(x,1)); 
         cellout{2}=adaptspec(abs(mmatx).^2,lambda,var).*dt;
     end
else         %Two time series
     if lambda==1
        cellout{2}=real(avgspec(mmatx,mmatx)).*dt;
        cellout{3}=real(avgspec(mmaty,mmaty)).*dt;
        cellout{4}=avgspec(mmatx,mmaty).*dt;
     else
        %For two time series one should do the adaptive spectra on both
        %with the same coefficients
        var=squared(vstd(x,1))+squared(vstd(y,1)); 
        
        [s,dk]=adaptspec(abs(mmatx).^2+abs(mmaty).^2,lambda,var);
        cellout{2}=squeeze(frac(sum(dk.^2.*abs(mmatx).^2.*dt,2),sum(abs(dk).^2,2)));
        cellout{3}=squeeze(frac(sum(dk.^2.*abs(mmaty).^2.*dt,2),sum(abs(dk).^2,2)));
        cellout{4}=squeeze(frac(sum(dk.^2.*mmatx.*conj(mmaty).*dt,2),sum(abs(dk).^2,2))); 
     end
end


%Corrections for zero component, and for Nyquist with even and odd length
%Both zero are Nyquist are shared for even length time series, so divide both by two.
%Only zero is shared for odd length, since the Nyquist does not appear. 
%This is best visualized by drawing N equally spaced points on the unit circle.

% for i=2:length(cellout)
%     cellout{i}(1,:)=cellout{i}(1,:)./2;
%     if iseven(size(x,1))
%        cellout{i}(end,:)=cellout{i}(end,:)./2;
%     end
% end

if strfind(normstr,'cyc')
    f=f/2/pi;
    %for i=2:length(cellout)
    %    cellout{i}=cellout{i}/2/pi;
    %end
end

    
    
if isscalar(dt)
    cellout{1}=f./dt;
else
    cellout{1}=vrep(f,size(dt,2),2)./dt;
end



function[S]=avgspec(mmat1,mmat2)
eigspec=mmat1.*conj(mmat2);
S=squeeze(mean(eigspec,2));



function[s,dk]=adaptspec(eigspec,lambda,var)

s=squeeze(zeros(size(eigspec(:,1,:))));
dk=zeros(size(eigspec));
for i=1:size(eigspec,3)
    [s(:,i),dk(:,:,i)]=adaptspec_one(eigspec(:,:,i),lambda,var(i));
end


function[s,dk]=adaptspec_one(eigspec,lambda,var)

s=frac(1,2)*squeeze(eigspec(:,1)+eigspec(:,2));
lambda=lambda(:)';
tol=1e-3;

var=var+0*s;
i=0;
sold=s*0+1000;
while any(abs(s-sold)./sold>tol)&&i<20
	i=i+1;
	sold=s;
    bk=var*(ones(size(lambda))-lambda);  %Outer product
	dk=(s*real(sqrt(lambda)))./(s*lambda+bk);  %Outer products
	s=frac(sum(dk.^2.*eigspec,2),sum(abs(dk).^2,2));
end

if i~=20
   disp(['Adaptive spectral estimate took ' int2str(i) ' iterations.'])
else
   disp(['Adaptive spectral loop terminated at ' int2str(i) ' iterations.'])
end


function[varargout]=mtrans(varargin)
% MTRANS  Multitaper "eigentransform" computation.
% 
%   [F,W]=MTRANS(X,PSI) returns the multitaper "eigentransform" matrix 
%   for use in multitaper spectral estimates or eigenspectral SVD analysis.                 
%
%       X  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
%       W  --  [M/2] x K x  N eigentransform matrix (for real X)
%              [M/2] x K x 2N eigentransform matrix (for complex X)
%
%   In the above, [M/2] means M/2 if M is even, and (M-1)/2 is M is odd.
%
%   F is the angular Fourier frequency, in radians.
%
%   [F,W1,W2,...,WN]=MTRANS(X1,X2,...,XN,PSI) also works, where X1,...,XN
%   are all M x N matrices.
%
%   MTRANS(X,[]) with PSI empty uses the default taper, so that the square
%   of W corresponds to the periodogram.  
%
%   See also: SLEPTAP, HERMFUN, MSPEC, MSVD.
%
%   Usage:  [f,w]=mtrans(x,psi);  
%           [f,wx,wy]=mtrans(x,y,psi);    
%           [f,wx,wy,wz]=mtrans(x,y,z,psi);    
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details        
  
%Sort out input arguments
psi=varargin{end};
x=varargin(1:end-1);

%Size check
for i=1:length(x)
    sizex(i,:)=size(x{i});
end
if ~allall(sizex(:,1)==sizex(1,1))||~allall(sizex(:,2)==sizex(1,2))
   error('All input arguments should be the same size')
end

sizex=sizex(1,:);

if ~isempty(psi)
    psimat=vrep(psi,sizex(2),3);
else
    psimat=[];
end

f=fourier(sizex(1));

index=1:length(f);
varargout{1}=f;

for i=1:size(x,2)
    mmat=mtrans1(x{i},psimat);
    varargout{i+1}=mmat(index,:,:);
end

function[mmat]=mtrans1(x,psimat)

Nnans=length(find(~isfinite(x)));
if Nnans>0
     disp(['MSPEC finding non-finite data values.  Spectrum will be undefined.'])
end

x=permute(x,[1 3 2]);
xmat=vrep(x,size(psimat,2),2);
mmat=fft(psimat.*xmat,[],1);

function[]=mspec_test

tol=1e-10;
load bravo94
cv=bravo.rcm.cv;
[psi,lambda]=sleptap(length(cv),8);

[f,sxx,syy,sxy]=mspec(real(cv),imag(cv),psi);
[f,spp,snn,spn]=mspec(cv,psi);

S(1,1,:,:)=sxx;
S(2,2,:,:)=syy;
S(1,2,:,:)=sxy;
S(2,1,:,:)=conj(sxy);

SZ(1,1,:,:)=spp;
SZ(2,2,:,:)=snn;
SZ(1,2,:,:)=spn;
SZ(2,1,:,:)=conj(spn);

T=vrep(vrep(tmat,size(S,3),3),size(S,4),4);

SZ2=matmult(matmult(T,S,1),conj(permute(T,[2 1 3 4])),1);

reporttest('MSPEC for (x,y) and (z,z^*) are unitary transforms with matrix T',aresame(SZ,SZ2,1e-10))


tol=1e-10;
x=bravo.rcm.cv(:,3);
xo=vfilt(x,24,'nonans');
%num=yf2num(bravo.rcm.yearf);
%t=num-yf2num(floor(bravo.rcm.yearf(1)));
t=bravo.rcm.yearf;


[psi,lambda]=sleptap(length(x),8);

p0=vsum(abs(real(x)-vmean(real(x),1)).^2,1)./length(x);
[f,sp]=mspec(real(x),psi);
p1=frac(1,2*pi)*(vsum(sp,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for real Bravo, unit sample rate',abs(p1-p0)./p0<4/100);

[f,sp]=mspec(t(2)-t(1),real(x),psi);
p1=frac(1,2*pi)*(vsum(sp,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for real Bravo, non-unit sample rate',abs(p1-p0)./p0<4/100);

[f,sp]=mspec(t(2)-t(1),real(x),[]);
p1=frac(1,2*pi)*(vsum(sp,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for real Bravo, non-unit sample rate periodogram',abs(p1-p0)./p0<4/100);

[f,sp]=mspec(t(2)-t(1),real(x),psi,'cyc');
p1=(vsum(sp,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for real Bravo, non-unit sample rate, cyclic frequency',abs(p1-p0)./p0<4/100);

p0=vsum(abs(x-vmean(x,1)).^2,1)./length(x);
[f,sp,sn,spn]=mspec(x,psi);
p1=frac(1,2*pi)*(vsum(sp,1)+vsum(sn,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for complex Bravo, unit sample rate',abs(p1-p0)./p0<4/100);

[f,sp,sn,spn]=mspec(t(2)-t(1),x,psi);
p1=frac(1,2*pi)*(vsum(sp,1)+vsum(sn,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for complex Bravo, non-unit sample rate',abs(p1-p0)./p0<4/100);

[f,sx,sy]=mspec(real(x),imag(x),psi);
p1=frac(1,2*pi)*(vsum(sx+sy,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for bivariate Bravo, unit sample rate',abs(p1-p0)./p0<4/100);

[f,sx,sy]=mspec(t(2)-t(1),real(x),imag(x),psi);
p1=frac(1,2*pi)*(vsum(sx+sy,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for bivariate Bravo, non-unit sample rate',abs(p1-p0)./p0<4/100);

reporttest('MSPEC SXX+SYY=SPP+SNN for Bravo, non-unit sample rate',aresame(sx+sy,sp+sn,tol));


N=1001;
cv=randn(N,1)+sqrt(-1)*randn(N,1);
cv=cv-vmean(cv,1);
[f,spp,snn]=mspec(cv,[]);
reporttest('MSPEC periodogram matches variance exactly for zero-mean complex signal and odd length', aresame(sum(spp)./N+sum(snn(2:end))./N,sum(abs(cv).^2)./N,1e-10))

N=1001;
cv=randn(N,1)+sqrt(-1)*randn(N,1)+17;
[f,spp,snn]=mspec(cv,[]);
reporttest('MSPEC periodogram matches variance exactly for non-zero-mean complex signal and odd length', aresame(sum(spp)./N+sum(snn(2:end))./N,sum(abs(cv).^2)./N,1e-10))


N=1000;
cv=randn(N,1)+sqrt(-1)*randn(N,1);
cv=cv-vmean(cv,1);
[f,spp,snn]=mspec(cv,[]);
reporttest('MSPEC periodogram matches variance exactly for zero-mean complex signal and even length', aresame(sum(spp)./N+sum(snn(2:end-1))./N,sum(abs(cv).^2)./N,1e-10))

N=1000;
cv=randn(N,1)+sqrt(-1)*randn(N,1)+17;
[f,spp,snn]=mspec(cv,[]);
reporttest('MSPEC periodogram matches variance exactly for non-zero-mean complex signal and even length', aresame(sum(spp)./N+sum(snn(2:end-1))./N,sum(abs(cv).^2)./N,1e-10))

N=1000;
dt=3600;
cv=randn(N,1)+sqrt(-1)*randn(N,1)+17;
[f,spp,snn]=mspec(dt,cv,[]);
reporttest('MSPEC periodogram matches variance exactly for non-zero-mean complex signal and even length, and non-unit sample rate', aresame(sum(spp)./(dt*N)+sum(snn(2:end-1))./(dt*N),sum(abs(cv).^2)./N,1e-10))

N=1000;
cv=randn(N,1)+sqrt(-1)*randn(N,1)+17;
tic;[f,spp,snn]=mspec(cv,[]);toc
tic;S=squared(fft(cv));toc;
b1=aresame(S(1:length(spp)),spp*N,1e-6);
b2=aresame(flipud(S(length(spp)+1:end)),snn(2:end-1)*N,1e-6);

reporttest('MSPEC recovers periodogram for even N', b1&&b2)

N=1001;
cv=randn(N,1)+sqrt(-1)*randn(N,1)+17;
tic;[f,spp,snn]=mspec(cv,[]);toc
tic;S=squared(fft(cv));toc;
b1=aresame(S(1:length(spp)),spp*N,1e-6);
b2=aresame(flipud(S(length(spp)+1:end)),snn(2:end)*N,1e-6);
reporttest('MSPEC recovers periodogram for odd N', b1&&b2)

mspec_test_frequency;

function[]=mspec_test_frequency


T=10000;
fo=100./T;
x=cos([1:T]'*2*pi*fo);
[psi,lambda]=sleptap(T,1,1);
[f,sp,sn]=mspec(x+sqrt(-1)./1e10,psi);

[mp,jp]=max(sp);
[mn,jn]=max(sn);
bool=aresame(frac(1,2*pi)*f(jp),fo,1e-12)&&aresame(frac(1,2*pi)*f(jn),fo,1e-12);
reporttest('MSPEC frequency matches expected exactly, even number of points',bool);

T=10000-1;
fo=100/T;
x=cos([1:T]'*2*pi*fo);
[psi,lambda]=sleptap(T,1,1);
[f,sp,sn]=mspec(x+sqrt(-1)./1e10,psi);

[mp,jp]=max(sp);
[mn,jn]=max(sn);
bool=aresame(frac(1,2*pi)*f(jp),fo,1e-12)&&aresame(frac(1,2*pi)*f(jn),fo,1e-12);
reporttest('MSPEC frequency matches expected exactly, odd number of points',bool);



% [f,Cuv]=mspec(real(cv),imag(x),psi);
% [f,Suu]=mspec(real(cv),psi);
% [f,Svv]=mspec(imag(cv),psi);
% 
% gammauv=Cuv;
% for i=1:size(Suu,2)
%   gammauv(:,i)=Cuv(:,i)./sqrt(Suu(:,i).*Svv(:,3));
% end
% figure,
% 
% 
% plot(f,abs(gammauv)),xlog,yoffset 1
% title('Cross-spectrum of u(t) at each depth vs. u(t) at #3')

