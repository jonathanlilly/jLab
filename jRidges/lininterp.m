function[x,t]=lininterp(t1,t2,x1,x2,t)
%LININTERP  Fast linear interpolation for arbitrary-sized arrays.
%
%   LININTERP is a low-level function called by RIDGEINTERP.
%
%   XI=LININTERP(T1,T2,X1,X2,TI) interpolates a function X within an 
%   independent variable T using a linear algorithm.
%
%   The function values X1 and X2 may be arrays of any size provided they
%   are the same size.  The "time" values T1 and T2 may be scalars, or 
%   arrays of the same size as X1 and X2.
%
%   If T1 and T2 and X1 and X2 are non-scalar arrays, T should either be
%   a scalar or an array of the same size as the other arguments.
%
%   If T1 and T2 and X1 and X2 are all scalars, T may be of any size, and
%   XI will have the same size as T.
%
%   The function values X1 and X2 may be real or complex.
%
%   LININTERP uses the exact algebraic expressions to operate looplessly
%   on matrices, thus it is very fast.    
%   __________________________________________________________________
%
%   Locating x-intercept
%
%   LININTERP can also be used to find the x-intercept value of X(T).
% 
%   TO=LININTERP(T1,T2,X1,X2) with only four input arguments returns TO,
%   the time at which X(T) vanishes.
%
%   TO is the same size as the input arrays.
%   __________________________________________________________________
%
%   Algorithm details
%
%   Let X and T be related as
%
%              X = A*T + B 
%
%   The coefficients A and B are uniquely determined by two (T,X) pairs.  
%   LININTERP solves for these coefficients and uses the result to 
%   interpolate X to any other value of T.   
%   __________________________________________________________________
%
%   See also QUADINTERP, CUBINTERP.
%  
%   Usage:   x=lininterp(t1,t2,x1,x2,t);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2012 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(t1, '--t')
    lininterp_test,return
end

if ~aresame(size(x1),size(x2))
    error('The input arrays X1 and X2 must be the same size.')
end

if ~aresame(size(t1),size(t2))
    error('The input arrays T1 and T2 must be the same size.')
end

maxlent=max([length(t1(:)) length(t2(:))]);
maxlenx=max([length(x1(:)) length(x2(:))]);

if nargin==5
    if maxlent>1||maxlenx>1
        if ~aresame(size(t1),size(t))&&~aresame(size(x1),size(t))&&~isscalar(t)
            error('The time T must either be the same size as the other input arguments, or a scalar.')
        end
    end
end

a=frac(x2-x1,t2-t1);
b=x1-a.*t1;

if nargin==4
    x=-frac(b,a); 
else
    x=a.*t+b;
end

function[]=lininterp_test

lininterp_test_real;
lininterp_test_complex;
lininterp_intercept_test_complex;

function[]=lininterp_test_real


a=1;
b=7;

t1=2;
t2=4;


%t=(-20:.1:20)';
t=[t1 t2];
x=a.*t+b;

N=100;
bool=false(N,1);
for i=1:100
   ti=randn(1)*10;
   xi=lininterp(t(1),t(2),x(1),x(2),ti);
   xi2=a.*ti+b;
   bool(i)=aresame(xi,xi2,1e-10);
end

reporttest('LININTERP real scalar case',all(bool(i)))


function[]=lininterp_test_complex
a=1+sqrt(-1)*4;
b=7+sqrt(-1)*2;

t1=2;
t2=4;

t=[t1 t2 ];
x=a.*t+b;

N=100;
bool=false(N,1);
for i=1:100
   ti=randn(1)*10;
   xi=lininterp(t(1),t(2),x(1),x(2),ti);
   xi2=a.*ti+b;
   bool(i)=aresame(xi,xi2,1e-10);
end

reporttest('LININTERP complex scalar case',all(bool(i)))


function[]=lininterp_intercept_test_complex
a=1+sqrt(-1)*4;
b=7+sqrt(-1)*2;

t1=2;
t2=4;

t=[t1 t2 ];
x=a.*t+b;

to=lininterp(t(1),t(2),x(1),x(2));
xo=lininterp(t(1),t(2),x(1),x(2),to);
 
reporttest('LININTERP x-intercept',aresame(xo,0))






