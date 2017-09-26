function[y,ystack]=slidetrans(x,w,fs,str)
%SLIDETRANS  Sliding-window ('moving-window') Fourier transform.
%   
%   Y=SLIDETRANS(X,PSI,F) computes the sliding-window Fourier transform of
%   the signal X using window PSI at frequencies F. X and W are column
%   vectors.  Y is a matrix of size LENGTH(X) by LENGTH(FS).
%
%   The frequencies F are radian frequencies as in COS(F T).   
%
%   PSI may also be a matrix, in which case Y is a 3-D array of size 
%   LENGTH(X) by LENGTH(FS) by SIZE(W,2).
%
%   SLIDETRANS(...,STR) sets the endpoint boundary conditions, with STR
%   equal to 'mirror', 'periodic', or 'zeros'.  The default is 'periodic'. 
%   For more details, see the boundary condition discussion in WAVETRANS.
%
%   A good choice for window is the first "Slepian" taper; see SLEPTAP.
%   Hermite functions are also used; see HERMFUN.
%
%   SLIDETRANS follows the same normalization as WAVETRANS when the input
%   signal X is complex-valued, so see that function for details.
%   _________________________________________________________________
%
%   Definition
%
%   The sliding window transform is defined as
%
%      y(t,f) = int psi^*(u-t) exp[-2 pi i f (u-t)] x(t) du.
%
%   When a sinusoid is transformed at its own frequency, the rate of change
%   of the phase of the sliding window transform recovers the frequency of
%   the sinusoid.  This is the same as with the wavelet transform.
%
%   This differs from the definition of Mallat (1999), p. 69, by a unit
%   amplitude factor of exp[2 pi i f t].  In Mallat's defintion, the
%   transform of a sinusoid at the sinusoid's frequency has constant phase.
%   Both definitions have the same modulus.
%   _________________________________________________________________
%
%   'slidetrans --t' runs a test.
%   'slidetrans --f' generates a sample figure.
%
%   Usage: y=slidetrans(x,psi,f);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details  

if strcmpi(x,'--f')
  type makefigs_slidetrans;
  makefigs_slidetrans
  return
end
if strcmpi(x,'--t')
  slidetrans_test;return
end

if nargin<4
    str='periodic';
end

M=size(x,1);
N=size(x,2);

fs=fs(:);
Mf=length(fs);
Mw=size(w,1);
K=size(w,2);

t=(0:M-1)';
tw=(0:Mw-1)';
tw=tw-mean(tw);


%Make a wavelet
w=vrep(permute(w,[1 3 2]),Mf,2);
phasor=vrep(rot(oprod(tw,fs)),K,3);
psi=w.*phasor;

%figure,uvplot(psi),hold on
%uvplot(psi,'o')

%psi=zeros(Mw,Mf,K);
%for j=1:Mf
%  for k=1:K
%    psi(:,j,k)=w(:,k).*rot(2*pi*tw*fs(j));
%  end
%end

y=wavetrans(x,psi,str);

%for j=1:Mf
%  for k=1:K
   % y(:,j,k)=y(:,j,k).*rot(-2*pi*t*fs(j)) ;
%  end
%end

if nargout==2
    test_freq(fs);
    ystack=0*y;
    for n=1:size(y,2)
        ystack(:,n,:)=vsum(y(:,n:n:end,:),2);
    end
end

function[]=test_freq(fs)

b=aresame(ones(length(fs)-1,1),diff(fs)./fs(1),1e-10);
if ~b
    error('Sorry, with "stack" option F must be of the form DF:DF:N*DF.')
end


function[]=slidetrans_test

N=501;
w=hermfun((-N:N)'./(N/4),0);
w=w./sum(w);
M=3001;
t=(0:M-1)';
x=rot(t./70);
y=slidetrans(x,w,1./70);
om=vdiff(unwrap(angle(y)),1);
bool(1)=aresame(om(N:M-N),1/70+0*(N:M-N)',1e-3);
bool(2)=aresame(abs(y(N:M-N)).^2,1/2+0*(N:M-N)',1e-1);
reporttest('SLIDETRANS complex sinusoid',all(bool))

% 
% 
% x=testseries(6);
% N=2000;
% w=hermfun([-N:N]'./(N/4),4);
% for i=1:size(w,2)
%   w(:,i)=w(:,i)./sqrt(w(:,i)'*w(:,i));
% end
% fs=(0:0.5:30]./length(x);
% y=slidetrans(x,w,fs);
% t=(0:length(x)-1)';
% jpcolor(t,fs,abs(y)'),shading interp,flipy
% jpcolor(t,fs,mean(abs(y(:,:,1:4)),3)'),shading interp,flipy
% contourf(t,fs,abs(y)',20),nocontours,flipy

function[]=slidetrans_fig_whale
  
%Not currently working
[x,fs,nbits]=audioread('callS1');
%[psi,lambda]=sleptap(size(x,1),16);
%[f,s]=mspec(x,psi);

x=x(1:10:end);

N=5000;
w=hermfun((-N:N)'./(N/4),0);
w=w./sqrt(w'*w);
fs=(100:100:length(w))./length(w);


[y,ystack]=slidetrans(x,w,fs,'zeros');
%y=slidetrans(x,w,fs,'zeros');

