function[varargout]=wavetrans(varargin)
%WAVETRANS  Continuous wavelet transform.
%   _______________________________________________________________________
%
%   *|* wavetrans.png --- Figure illustrating the analytic wavelet 
%   transform. Type 'jhelp wavetrans' to view this image. *|*
%   _______________________________________________________________________
%
%   W=WAVETRANS(X,PSI) computes the wavelet transform W of a data series X
%   using wavelets PSI. X has time oriented in rows, and PSI is either a
%   matrix of wavelets, or a cell array defining wavelet properties.
%
%   If PSI is a matrix containing the wavelets, SIZE(PSI,1) should be the 
%   same as SIZE(X,1), with time-domain wavelets at different scales or 
%   frequencies in the columns of PSI.  See MORSEWAVE or MORLWAVE.
%
%   If PSI is a cell array, it defines the properties of generalized Morse
%   wavelets, which are then computed internally, as described below. 
%   ___________________________________________________________________
%
%   Multiple input series
%
%   X and PSI may both contain multiple components.  Then PSIF(:,:,k) is
%   the kth wavelet, and X(:,n) is the nth data component.  If there are K
%   wavelets at J frequencies and M data points in N components of X, then
%   W is of size M x J x N x K.  Note that W is always squeezed to remove
%   singleton dimensions.
%
%   X can also be a 3D array of size M x N x K, if the wavelet contains
%   only one component.  Then W is again of size M x J x N x K.  
%
%   [W1,W2,...,WN]=WAVETRANS(X1,X2,...,XN,PSI) also works, where the XN are
%   all column vectors of the same length.
%   ___________________________________________________________________
%
%   Generalized Morse wavelets
%
%   WAVETRANS can automatically compute the wavelet transform using the
%   generalized Morse wavelets, without needing to precompute the wavelets.
%
%   WAVETRANS(X,{GAMMA,BETA,FS}), with PSI being a cell array, uses the
%   Generalized Morse Wavelets specified by the parameters GAMMA and BETA.
%
%   FS is an array of *radian* frequencies, as in cos(ft) not cos(2 pi ft),
%   assuming a unit sample rate.  Thus the Nyquist frequency is at pi.
%   Use MORSESPACE to easily choose the frequency bins.
%
%   WAVETRANS(X,{K,GAMMA,BETA,FS}) uses the first K orthogonal multi-
%   wavelets to create K different transforms.  The default is K=1,
%   employing only the first wavelet.  See MORSEWAVE for details.
%
%   WAVETRANS(X,{...,'energy'}) specifies a unit energy normalization for
%   the wavelets.  The default is the 'bandpass' or unit amplitude
%   normalization.  See MORSEWAVE for details.
%
%   For general purpose use, set GAMMA=3 and choose BETA to be no smaller
%   than one.  Increase BETA to make your wavelet have more 'wiggles'.
%   ___________________________________________________________________
%
%   Boundary conditions
%
%   W=WAVETRANS(...,STR), where STR is a string, optionally specifies the
%   boundary condition to be imposed at the edges of the time series.
%   Valid options for STR are
%
%         STR = 'periodic' for periodic boundary conditions
%         STR = 'zeros' for zero-padding beyond the endpoints
%         STR = 'mirror' for reflecting the time series at both ends
%         STR = 'reverse' for reflection together with a sign reversal
%
%   The default value of STR is 'periodic', which means endpoints of the
%   time series are implicitly joined to make a periodic signal. All
%   boundary conditions take into account potential blocks of missing data,
%   marked by NaNs, at beginning and end of each column.
%   ___________________________________________________________________
%
%   Missing data
%
%   The data X may contain blocks of NANs at the beginning and/or end of
%   each column, marking the absence of data.  In this case only the
%   data series is taken to correspond to the block of finite data values,
%   and the boundary conditions are applied accordingly. The corresponding
%   portions of the transform matrix W are then also set to NANs. No NANs
%   may occur in the interior of the data series.
%   ___________________________________________________________________
%
%   Detrending
%
%   Note that the data X is detrended before transforming.  This feature
%   is suppressed by WAVETRANS(..., 'nodetrend').
%   ___________________________________________________________________
%
%   Complex-valued data
%
%   The wavelet transform is normalized differently for complex-valued data
%   than for real-valued data, and this in turns depends on whether the 
%   'bandpass' or 'energy' normalizations are specified. 
%
%   If WX and WY are the wavelet transforms of two real-valued signals, 
%   X and Y, then if the 'bandpass' normalization is specified, 
%
%        WP=WAVETRANS(X+iY,PSI)   = (1/2)*(WX + i WY)
%        WN=WAVETRANS(X-iY,PSI)   = (1/2)*(WX - i WY)
%
%   defines the positive and negative rotary transforms WP and WN.
%
%   The factor of 1/2 sets the peak value of the wavelet transform of a
%   complex exponential to the amplitude of the complex exponential itself.
%
%   Alternatively, if the 'energy' normalization is specified, 
%
%        WP=WAVETRANS(X+iY,PSI)   = (1/SQRT(2))*(WX + i WY)
%        WN=WAVETRANS(X-iY,PSI)   = (1/SQRT(2))*(WX - i WY)
%
%   is used instead. The factors of SQRT(2) are included such that the
%   total power is unchanged between the rotary and Cartesian transforms:
%
%        ABS(WX).^2+ABS(WY).^2 = ABS(WP).^2+ABS(WN).^2.
%
%   There are two equivalent ways to compute the positive and
%   negative rotary transforms.  The first is
%
%     [WP,WN]=WAVETRANS(X+iY,X-iY,PSI) 
%
%   which takes the wavelet transform of the X+iY and its conjugate.
% 
%   The second is to take the wavelet transform of the real-valued data
%
%       [WX,WY]=WAVETRANS(X,Y,PSI)              followed by 
%       [WP,WN]=VECTMULT(TMAT/SQRT(2),WX,WY)    (bandpass normalization) or 
%       [WP,WN]=VECTMULT(TMAT,WX,WY)            (energy normalization) 
%
%   which converts WX and WY to WP and WN with a matrix multiplication
%   using the unitary matrix TMAT=[1 i; 1 -i]/SQRT(2).
%   ___________________________________________________________________
%
%   Parallelization
%
%   WAVETRANS(...,'parallel') parallelizes the wavelet transform using a 
%   PARFOR loop, when the data X contains multiple columns.  This requires 
%   that Matlab's Parallel Computing Toolbox be installed. 
%   __________________________________________________________________
%
%   See also MORSESPACE, RIDGEWALK, WAVESPECPLOT.
%
%   'wavetrans --t' runs some tests.
%   'wavetrans --f' generates a sample figure.
%
%   Usage:  w=wavetrans(x,psi);
%           w=wavetrans(x,{gamma,beta,f,str});
%           w=wavetrans(x,{gamma,beta,f,str},str);
%           [wx,wy]=wavetrans(x,y,{gamma,beta,f,str},str);
%           [wp,wn]=wavetrans(x+i*y,x-i*y,{gamma,beta,f,str},str);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2021 J.M. Lilly --- type 'help jlab_license' for details


if strcmpi(varargin{1},'--t')
    wavetrans_test;return
elseif strcmpi(varargin{1},'--f')
    type makefigs_wavetrans 
    makefigs_wavetrans;
    return
end

%Look for 'parallel' or 'series' input 
parstr='serial';
for i=1:length(varargin)
    if ischar(varargin{i})
        if strcmpi(varargin{i}(1:3),'ser')||strcmpi(varargin{i}(1:3),'par')
            parstr=varargin{i};
            varargin=varargin([1:i-1 i+1:length(varargin)]);
        end
    end
end

bool=zeros(length(varargin),1);
for i=1:length(varargin)
    bool(i)=ischar(varargin{i});
end
if ~allall(bool==0)
    N=find(bool,1,'first')-2;
else
    N=length(varargin)-1;
end
argcell=varargin(N+1:end);
for j=1:N
    x=varargin{j};
    if iscell(x)
        if strcmpi(parstr(1:3),'par')
            disp('WAVETRANS employing parallel algorithm.')
            parfor i=1:length(x)
                disp(['WAVETRANS transforming time series ' int2str(i) ' of ' int2str(length(x)) '.'])
                %[T{i},fs{i}]=wavetrans_one(x{i},argcell);
                T{i}=wavetrans_one(x{i},argcell,'serial');
            end
        else
            for i=1:length(x)
                disp(['WAVETRANS transforming time series ' int2str(i) ' of ' int2str(length(x)) '.'])
                %[T{i},fs{i}]=wavetrans_one(x{i},argcell);
                T{i}=wavetrans_one(x{i},argcell,'serial');
            end
        end
    else
        %[T,fs]=wavetrans_one(x,argcell);
        T=wavetrans_one(x,argcell,parstr);
    end
    varargout{j}=T;
end

% if ~isempty(fs)
%     varargout{j+1}=fs;
% end

function[T]=wavetrans_one(x,argcell,parstr)
%This is just so we don't bother with empty columns

xpages=size(x,3);
if xpages>1
    x=reshape(x,[size(x,1),size(x,2)*size(x,3)]);
end

ngood=sum(isfinite(x)+0,1);
goodindex=find(ngood>1);

if isempty(goodindex)
    %Just pass one through so we can get the right size
    T=wavetrans_continue(zeros(size(x(:,1))),argcell,parstr);
    T=inf*(1+sqrt(-1)).*vrep(T,size(x,2),2);
else
    x=x(:,goodindex);
    T1=wavetrans_continue(x,argcell,parstr);
    T=inf*(1+sqrt(-1)).*zeros(size(T1,1),size(T1,2),size(x,2),size(T1,4));
    T(:,:,goodindex,:)=T1;
end

if xpages>1
    T=reshape(T,[size(T,1) size(T,2) size(T,3)./xpages xpages]);
end

T=squeeze(T);


function[T]=wavetrans_continue(x,argcell,parstr)

normstr='bandpass'; %default for wavetrans
if iscell(argcell{1})
    if ischar(argcell{1}(end))
        normstr=argcell{1}(end);
    end
end

if ~isreal(x)
    %Unitary transform normalization --- use this only for 'energy'
    if strcmpi(normstr(1:3),'ene')
        x=x./sqrt(2);
    %Bandpass transform normalization --- use this for 'bandpass'
    elseif strcmpi(normstr(1:3),'ban')
        x=x./2;
%        x=x./sqrt(2); %XXXX
    end
end

w=argcell{1};
detrendstr='detrend';
str='periodic';

if ischar(argcell{end})&&~isempty(argcell{end})
    if strcmpi(argcell{end}(1:3),'nod')
        detrendstr='nodetrend';
        if ischar(argcell{end-1})
            str=argcell{end-1};
        end
    else
        str=argcell{end};
    end
end


x0=x;
M0=size(x0,1);

if isreal(x0)
    %figure,plot(x0)
    x=timeseries_boundary(x0,str,detrendstr);
    %hold on,plot(x,'r')
else
    x=           timeseries_boundary(real(x0),str,detrendstr)+...
        sqrt(-1)*timeseries_boundary(imag(x0),str,detrendstr);
end
M=size(x,1);
N=size(x,2);
W=[];

%figure,plot(x)
if ~iscell(w)
    if size(w,1)>M0 && strcmpi(str,'periodic')
        disp('Data length must exceed the filter length---returning INFs.')
    elseif size(w,1)>3*M0
        disp('Data length must exceed one-third the filter length---returning INFs.')
    end
end


%/********************************************************
%fs=[];
if iscell(w)
    %     if length(w)==3
    %         if ischar(w{end})
    %             w{4}=w{end};
    %             w{3}=morsespace(w{1},w{2},M);
    %         end
    %     elseif length(w)==2
    %         w{3}=morsespace(w{1},w{2},M);
    %     end
    [psi,W]=morsewave(M,w);
    K=size(W,3);
    L=size(W,2);
    %    fs=w{3};
else
    K=size(w,3);
    L=size(w,2);
    
    %Generate a frequency-domain wavelet matrix of same size as data
    if size(w,1)<M   %Some subtlety here for even/odd or odd/even
        wnew=zeros(M,L,K);
        index=(1:size(w,1))'+floor((M-size(w,1))./2);
        wnew(index,:,:)=w;
        w=wnew;
    elseif size(w,1)>M
        w=nan*w(1:M,:,:,:);  %If wavelet is too long, just truncate
    end
    
    W=fft(w);
    om=2*pi*linspace(0,1-1./M,M)';
    om=vrep(om,L,2);
    om=vrep(om,K,3);
    
    if iseven(M)
        W=W.*rot(-om.*(M+1)/2).*sign(pi-om);  %ensures wavelets are centered
    else
        W=W.*rot(-om.*(M+1)/2);               %ensures wavelets are centered
    end
    %Note, the sign function you need when the wavelets are real-valued for
    %even-length time series only
    
    W=reshape(W,M,L,K);
end
%\********************************************************

%figure,plot(ifft(W))
%figure,plot(imag(W))
W=conj(W);
%figure,plot(abs(W))
X=fft(x);
T=nan*ones(M0,L,N,K);
%figure,plot(abs(W))

if M0~=M
    index=M0+1:M0*2;
else
    index=(1:M0);
end


if strcmpi(parstr(1:3),'par')
    disp('WAVETRANS employing parallel algorithm.')
    for k=1:K
        parfor n=1:N;
            Ttemp=ifft(vrep(X(:,n),L,2).*W(:,:,k));
            T(:,:,n,k)=Ttemp(index,:);
        end
    end
else
    for k=1:K
        for n=1:N;
            Ttemp=ifft(vrep(X(:,n),L,2).*W(:,:,k));
            T(:,:,n,k)=Ttemp(index,:);
        end
    end
end


%figure,plot(abs(fft(T)),'b')
%figure,plot(T,'m')

if isreal(x0)&&~iscell(w)
    if isreal(w)
        if ~isreal(T)
            T=real(T);  %Strip imaginary part if real wavelet and real signal
        end
    end
end

%figure,plot(abs(fft(T)),'m')
if ~anyany(isfinite(T))
    if ~isreal(T)
        T=inf*(1+sqrt(-1))*ones(size(T));
    else
        T=inf*ones(size(T));
    end
end

%/********************************************************
%Set missing data back to NANs
for i=1:size(x,2)
    index=find(isnan(x0(:,i)));
    if~isempty(index)
        Ttemp=T(:,:,i);
        if ~isreal(T)
            Ttemp(index,:)=nan*(1+sqrt(-1));
        else
            Ttemp(index,:)=nan;
        end
        
        T(:,:,i)=Ttemp;
    end
end
%\********************************************************


function[]=wavetrans_test
wavetrans_test_centered;
wavetrans_test_sizes;
wavetrans_test_complex;
wavetrans_test_boundary;
wavetrans_test_tooshort;

function[]=wavetrans_test_tooshort
load npg2006
%use npg2006
cx=npg2006.cx;

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
psi=morsewave(length(cx)+10,1,2,4,fs,'bandpass');
wx=wavetrans(real(cx),psi,'periodic');
reporttest('WAVETRANS returns INFs when signal is too short, periodic',allall(~isfinite(wx)))
psi=morsewave(3*length(cx)+10,1,2,4,fs,'bandpass');
wx=wavetrans(real(cx),psi,'mirror');
reporttest('WAVETRANS returns INFs when signal is too short, mirror',allall(~isfinite(wx)))

% function[]=wavetrans_test_energy
%
% fs=2*pi./(logspace(log10(10),log10(100),50)');
% cx=morsewave(1024,1,2,4,fs(round(end/2)),'energy');
%
% wx=wavetrans(frac(2,sqrt(2))*real(cx),{1,2,4,fs,'energy'},'mirror');
% reporttest('WAVETRANS unit energy',aresame(maxmax(abs(wx)),1,1/20))


function[]=wavetrans_test_boundary
load npg2006
%use npg2006
cx=npg2006.cx;
cv=npg2006.cv;

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
wx=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'mirror');
wx2=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'periodic');
wx3=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'zeros');

res1=abs(wx-wx2);
res2=abs(wx-wx3);

reporttest('WAVETRANS mirror and periodic boundary conditions match in interior',allall(res1(200:920)<1e-2))
reporttest('WAVETRANS mirror and zero boundary conditions match in interior',allall(res2(200:920)<1e-2))

% cx=cx(1:end-3);
% wx=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'mirror');
% [psi,psif]=morsewave(length(cx),1,2,4,fs);
% wx4=wavetrans(real(cx),psi,'mirror');
% reporttest('WAVETRANS input wavelet matches with automatic computation',aresame(wx,wx4,1e-10))

wx=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'periodic');
[psi,psif]=morsewave(length(cv),1,2,4,fs);
wx4=wavetrans(real(cx),psi,'periodic');
reporttest('WAVETRANS input wavelet matches with automatic computation',aresame(wx,wx4,1e-10))
wx5=wavetrans(real(cx),real(psi),'periodic');
reporttest('WAVETRANS real part of transform equals transform with real wavelet, odd length',aresame(real(wx4),wx5,1e-10))

vindex(cx,cv,1:1116,1);
[psi,psif]=morsewave(length(cv),1,2,4,fs);
wx4=wavetrans(real(cx),psi,'periodic');
wx5=wavetrans(real(cx),real(psi),'periodic');
reporttest('WAVETRANS real part of transform equals transform with real wavelet, even length',aresame(real(wx4),wx5,1e-10))



function[]=wavetrans_test_complex
load npg2006
%use npg2006
cx=npg2006.cx;
%cv=npg2006.cv;

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp,wn]=wavetrans(cx,conj(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp2,wn2]=vectmult(tmat./sqrt(2),wx,wy);


reporttest('WAVETRANS complex-valued input',aresame(wp,wp2,1e-6)&&aresame(wn,wn2,1e-6))

function[]=wavetrans_test_sizes

x=testseries_lillypark_array;
M=size(x,1);
N=size(x,2);
%Calculate wavelet matrix
J=5;
K=3;
fs=1./(logspace(log10(20),log10(600),J)');
psi=morsewave(M,K,2,4,fs,'bandpass');
%Compute wavelet transforms
wx=wavetrans(x,psi,'mirror');
[M2,J2,N2,K2]=size(wx);

bool=aresame([M,J,K,N],[M2,J2,K2,N2]);

reporttest('WAVETRANS output matrix has size M x J x N x K',bool)

function[]=wavetrans_test_centered


J=4;
ao=logspace(log10(5),log10(40),J)/100;
x=zeros(2^10-1,1);t=(1:length(x))';
[w,f]=morsewave(length(x),1,2,4,ao);
x(2^9,1)=1;
y=wavetrans(x,w);
clear maxi
for i=1:size(y,2);
    maxi(i)=find(abs(y(:,i))==max(abs(y(:,i))));
end
b(1)=max(abs(maxi-2^9)<=1);
reporttest('WAVETRANS Morsewave transform has peak at delta-function',b(1))


function[x,t]=testseries_lillypark_array
[x1,t]=testseries_lillypark;
x1=anatrans(x1);
x1=x1(1:10:end);
t=t(1:10:end);
dom=pi/6;
p1=[1 2 3 3 2.5 1.5]'.*rot(dom*(0:5)');
p1=p1./sqrt(p1'*p1);
x=real(oprod(x1,p1));
x=10*x./maxmax(x);

function[x,t,xo]=testseries_lillypark
t=(1:4096)';
om=zeros(size(t));
x2=0*om;
index=1000:3000;
om(index)=(t(index))*2/100000 - .019999;
x=sin(om.*t);
x2(index)=-sin((t(index)-1600)*2*pi*1/2800);
x=x.*x2;
x(3500:3515)=-1;
x(3516:3530)=1;
t=linspace(-100,100,length(x))';
xo=x2;
