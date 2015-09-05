function[varargout]=blurspec(varargin)
%BLURSPEC  Returns the blurred and aliased spectrum given the autocovariance.
%
%   BLURSPEC is used to rapidly compute blurred, aliased, and other 
%   modified versions of a spectrum from a known autocovariance.
%
%   Performing these calculations in the time domain, making use of 
%   Fourier relationships between modified versions of the autocovariance
%   and spectrum, is much faster than working in the frequency domain.
%   __________________________________________________________________
%
%   Blurred and aliased spectra
%
%   [F,S]=BLURSPEC(R) inverse Fourier transforms a given *one-sided*
%   autocovariance R to obtain a spectrum S that incorporates the effects
%   of aliasing, as well as blurring by the default taper, the 'boxcar'.
%
%   Here F is an array of frequencies and S is the one-sided Fourier 
%   spectrum.  If R is length N, the output arrays will be length (N/2+1) 
%   if N is even and (N+1)/2 if N is odd.
%
%   [F,SPP,SNN]=BLURSPEC(R) for complex-valued R returns SPP and SNN, the 
%   positive and negative rotary spectra, respectively.  
%
%   BLURSPEC(DT,R) optionally uses the sample interval DT in computing the 
%   frequency array F, and in setting the spectral levels.
%   __________________________________________________________________
%
%   Tapered and aliased spectra
%
%   BLURSPEC(R,'tapered',TAPER) incorporates the spectral smoothing due to 
%   the use of data TAPER, rather than the boxcar taper, in addition to the
%   effects of aliasing.  TAPER must be the same length as R.
%
%   Note that BLURSPEC(R,'tapered',[]) simply returns the blurred spectrum,
%   as an empty taper is taken to indicate the periodogram.
%
%   BLURSPEC(R,'window',WIN) uses the pre-computed window that is to 
%   multiply the autocovariance function WIN.  This is the half of the 
%   sequence obtained by convolving the taper with itself.  This version
%   is primarily used for speed in an internal call from MATERNFIT.
%   __________________________________________________________________
%
%   Aliased-only spectrum
%
%   BLURSPEC(R,'aliased') computes an approximation to *aliased* spectrum,
%   without blurring.  This approximation will be accurate to the extent
%   that R has decayed to zero by the end of its duration.
%
%   This is much faster than explicitly summing over aliased frequencies.
%   __________________________________________________________________
%
%   Differenced spectra
%
%   BLURSPEC(R,'difference') returns the blurred and aliased spectrum of 
%   the first forward difference of the process with autocovariance R.  In 
%   this case the output will be length N-1.
%
%   BLURSPEC(R,'seconddifference') returns the blurred and aliases spectrum
%   of the second forward difference of the process with autocovariance R.  
%   The output will be length N-2.
%
%   These options can be combined with the 'taper' and 'aliased' options 
%   described above. 
%   __________________________________________________________________
%
%   'blurspec --t' runs some tests.
%
%   Usage: [f,S]=blurspec(R);
%          [f,Spp,Snn]=blurspec(dt,R);        
%          [f,Spp,Snn]=blurspec(dt,R,'tapered',TAPER);        
%          [f,Spp,Snn]=blurspec(dt,R,'aliased');        
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly and A.M. Sykulski 
%                                 --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    blurspec_test,return
end

dt=1;
if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
end
R=varargin{1};

ver='standard';
str='blurred';
psi=[];

varargin=varargin(2:end);
for i=1:2
    if length(varargin)>0
        if ischar(varargin{end})
            if strcmpi(varargin{end}(1:3),'sec')||strcmpi(varargin{end}(1:3),'dif')||strcmpi(varargin{end}(1:3),'sta')
                ver=varargin{end};
            elseif strcmpi(varargin{end}(1:3),'ali')
                str=varargin{end};
            end
            varargin=varargin(1:end-1);
        end
    end
    if length(varargin)>1
        %length(varargin)
        if ischar(varargin{end-1})
            str=lower(varargin{end-1});
            psi=varargin{end};
            varargin=varargin(1:end-2);
        end
    end
end

if isempty(psi)&&strcmpi(str(1:3),'tap')
    str='blurred';
end

if strcmpi(ver(1:3),'dif')
    R=frac(1,dt)*(2*R(1:end-1,:)-R(2:end,:)-[conj(R(2,:));R(1:end-2,:)]);
elseif strcmpi(ver(1:3),'sec')
    R=frac(1,dt)*(2*R(1:end-1,:)-R(2:end,:)-[conj(R(2,:));R(1:end-2,:)]);
    R=frac(1,dt)*(2*R(1:end-1,:)-R(2:end,:)-[conj(R(2,:));R(1:end-2,:)]);
end

N=size(R,1);
if strcmpi(str(1:3),'win')
    R=R.*vrep(psi,size(R,2),2);
elseif strcmpi(str(1:3),'tap')
    if length(psi)~=size(R,1)
        error('Sizes of R and taper PSI do not match.')
    end
    win=conv(psi,psi);
    win=win(end-N+1:end);
    R=R.*vrep(win,size(R,2),2);
elseif strcmpi(str(1:3),'blu')
    tri=[N:-1:1]'./N;
    if size(R,2)==1
        R=R.*tri;
    else
        R=R.*vrep(tri,size(R,2),2);
    end
end  %If str='aliased', do nothing
%figure,plot(tri)

R(1,:)=R(1,:)./2;  %Don't forget to divide first element by two
S=dt*2*real(fft(R));
S=abs(S);  %Sometimes there are small negative parts after blurring

%Note, this is always correct for both even and odd length time series
omega=fourier(N);
Spp=S(1:length(omega),:);
Snn=[S(1,:);S(end:-1:end-length(omega)+2,:)];
%Snn=flipud([S(end-length(omega)+2:end,:);S(1,:)]);  %Same but slower

varargout{1}=omega./dt;
varargout{2}=Spp;
varargout{3}=Snn;

function[]=blurspec_test
 
N=1000;
alpha=1.5;
h=1;
A=1;

[tau,R]=materncov(1,N,A,alpha,h);
[f,Spp,Snn]=maternspec(1,N,A,alpha,h);
tic;[f,Spp2,Snn2]=blurspec(R,'aliased');etime1=toc;

maternspec_spec=@(omega,A,omegao,H,alpha)(frac(H.^(2*alpha-1),materncfun(alpha))*frac(A.^2,((omega-omegao).^2+H.^2).^alpha));
tic;
M=100;
[Sppa,Snna]=vzeros(length(fourier(N)),2*M+1);
for m=-M:M
    Sppa(:,m+M+1)=maternspec_spec(fourier(N)+2*pi*m,A,0,h,alpha);
end
Snna=Sppa;
etime2=toc;

vsum(Sppa,Snna,2);
%figure,plot(f,[Spp2 Spp2],'k','linewidth',2),hold on,plot(f,[Sppa Snna],'r')
%figure,plot(f,Spp2-Sppa,'k','linewidth',2),hold,plot(f,Snn2-Snna,'r','linewidth',2)

b1=aresame(Spp2,Sppa,2e-6);
b2=aresame(Snn2,Snna,2e-6);

reporttest('BLURSPEC aliasing with covariance inversion matches direct calculation to 2e-6, even N',b1&&b2)
disp(['BLURSPEC aliasing with covariance inversion was ' num2str(etime2/etime1) ' times faster than direct calculation.'])


N=999;
omegao=1/4;
[tau,R]=materncov(1,N,A,alpha,h,omegao);
[f,Spp,Snn]=maternspec(1,N,A,alpha,h,omegao);
tic;[f,Spp2,Snn2]=blurspec(R,'aliased');etime1=toc;

tic;
M=200;
[Sppa,Snna]=vzeros(length(fourier(N)),2*M+1);
for m=-M:M
    Sppa(:,m+M+1)=maternspec_spec(fourier(N)+2*pi*m,A,omegao,h,alpha);
    Snna(:,m+M+1)=maternspec_spec(fourier(N)+2*pi*m,A,-omegao,h,alpha);
end
etime2=toc;

vsum(Sppa,Snna,2);
b1=aresame(Spp2,Sppa,1e-6);
b2=aresame(Snn2,Snna,1e-6);
  
reporttest('BLURSPEC aliasing with covariance inversion matches direct calculation to 1e-6, odd N, frequency shift',b1&&b2)
disp(['BLURSPEC aliasing with covariance inversion was ' num2str(etime2/etime1) ' times faster than direct calculation.'])

x=[10 1.1 0.1];
[tau,acv]=materncov(1,N,x(1),x(2),x(3)); % autocovariance sequence
S=abs(real(2*fft(acv.*(1-([0:N-1]')/N))-acv(1))); % blurred spectrum

Spp2=S(1:floor(N/2)+1);
Snn2=[S(1);S(end:-1:ceil(N/2)+1)]; 
[tau,R]=materncov(1,N,x(1),x(2),x(3));
[f,Spp,Snn]=blurspec(R);

%[f,Spp,Snn]=maternspec(N,x(1),x(2),x(3),'blurred');

reporttest('BLURSPEC blurring matches alternate version for Matern', aresame(Spp,Spp2,1e-8)&&aresame(Snn,Snn2,1e-8))

% N=1000;
% alpha=10;
% h=1/10;
% sigma=7;
% z=maternoise(1,[N 10000],sigma,-1/2,h,0,alpha);
% psi=sleptap(size(z,1),3,1);
% [f,spp,snn]=mspec(1,z,psi);
% 
% 
% [tau,R]=materncov(1,N,sigma,-1/2,h,0,alpha);
% [f,Spp,Snn]=blurspec(R,'taper',psi);
% figure,plot(f,vmean([spp snn],2)),hold on,plot(f,Spp)
% 
% N=1001;
% z=maternoise(1,[N 10000],sigma,-1/2,h,0,alpha);;
% psi=sleptap(size(z,1),3,1);
% [f,spp,snn]=mspec(1,z,psi);
% 
% 
% [tau,R]=materncov(1,N,sigma,-1/2,h,0,alpha);
% [f,Spp,Snn]=blurspec(R,'taper',psi);
% %figure,plot(f,vmean([spp snn],2)),hold on,plot(f,Spp)
% 
% 
% N=1000;
% alpha=1.5;
% h=1;
% A=1;
% psi=sleptap(size(z,1),3,1);
% 
% [f,Spp,Snn]=maternspec(dt,N,A,alpha,h/dt);
% [tau,R]=materncov(dt,N,A,alpha,h);
% tic;[f,Spp2,Snn2]=blurspec(dt,R,'tapered',psi);
% tic;[f,Spp3,Snn3]=blurspec(dt,R);
% figure,plot(f,[Spp Spp2 Spp3])
% 
% 
% dt=7;
% z=maternoise(dt,[N 1000],A,alpha,h/dt);
% psi=sleptap(size(z,1),3,1);
% [f,spphat,snnhat,spn]=mspec(dt,z,psi);
% vmean(spphat,snnhat,2);
% 
% 


