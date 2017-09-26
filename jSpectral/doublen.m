function[y]=doublen(x)
% DOUBLEN  Interpolates a time series to double its length.
%
%   Y=DOUBLEN(X) for a column vector X of length M, returns the
%   interpolation of X to length 2M.  M must be even.
%  
%   For matrix X of size M x N, each column is interpolated to length
%   2M, and Y is size 2M x N.
%  
%   The interpolated vector Y may be considered to be of the same
%   duration, but with twice the sample rate, as the original vector.
%   In the interpolated vector, all frequencies higher than the
%   original Nyquist are set to zero.
%  
%   DOUBLEN uses the frequency-domain zero padding algorithm of Mallat
%   (1999), second edition, Section 4.5.4.
%
%   'doublen --t' runs a test.
%   'doublen --f' generates a sample figure.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015  J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(x,'--t')
    doublen_test;return
end

if strcmpi(x,'--f')
    type makefigs_doublen
    makefigs_doublen;return
end
 
M=size(x,1);

if isodd(M)
  error('The number of rows of X must be even.')
end
N=size(x,2);
xhat=fft(x);
yhat=zeros(2*M,N);

%Positive frequencies
index=(1:M/2);
yhat(index,:)=2*xhat(index,:);

%Negative frequencies
index=(3*M/2:2*M);
yhat(index,:)=2*xhat(index-M,:);

%Nyquist (twice)
yhat(M/2,:)=xhat(M/2,:);
yhat(3*M/2,:)=xhat(M/2,:);

y=ifft(yhat);
if allall(isreal(x))
   %strip small imaginary part for real x
   y=real(y);
end
   

function[]=doublen_test

M=100;
N=100;
x=randn(M,N);
y=doublen(x);

t1=(0:100-1)';
t2=(0:1/2:100-1/2)';

tol=1e-10;
b=aresame(y(1:2:end),x(1:end),tol);
strM= int2str(M);
strN= int2str(N);

reporttest(['DOUBLEN Y(1:2:end,:)=X for ' strM  ' x ' strN ' random trial'],b);

