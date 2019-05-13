function[x,t]=quadinterp(t1,t2,t3,x1,x2,x3,t)
%QUADINTERP  Fast quadratic interpolation for arbitrary-sized arrays.
%
%   QUADINTERP is a low-level function called by RIDGEINTERP.
%
%   XI=QUADINTERP(T1,T2,T3,X1,X2,X3,TI) interpolates a function X
%   within an independent variable T using a quadratic algorithm.
%
%   The function values X1--X3 may be arrays of any size provided 
%   they are all the same size.  The "time" values T1--T3 may be 
%   scalars, or arrays of the same size as X1--X3.  
%
%   If T1--T3 and X1--X3 are non-scalar arrays, T should either be
%   a scalar or an array of the same size as the other arguments.
%
%   If T1--T3 and X1--X3 are all scalars, T may be of any size, and
%   XI will have the same size as T.
%
%   The function values X1--X3 may be real or complex.
%
%   QUADINTERP uses the exact algebraic expressions to operate
%   looplessly on matrices, thus it is very fast.    
%   __________________________________________________________________
%
%   Locating maxima or minima
%
%   QUADINTERP can also be used to find the extreme value of X(T), 
%   that is its maxima or minima.
% 
%   [XE,TE]=QUADINTERP(T1,T2,T3,X1,X2,X3) with only six input arguments
%   solves for TE, the time at which X(T) obtains an extreme value XE.
%
%   XE and TE are the same size as the input arrays.
%   __________________________________________________________________
%
%   Algorithm details
%
%   Let X and T be related as
%
%              X = A*T.^2 + B*T + C 
%
%   The coefficients A, B, and C are uniquely determined by three (T,X)
%   pairs.  QUADINTERP solves for these coefficients and uses the 
%   result to interpolate X to any other value of T.   
%   __________________________________________________________________
%
%   See also LININTERP, CUBEINTERP.  
%  
%   Usage:   x=quadinterp(t1,t2,t3,x1,x2,x3,t);
%            [xe,te]=quadinterp(t1,t2,t3,x1,x2,x3,t);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(t1, '--t')
    quadinterp_test,return
end

if ~aresame(size(x1),size(x2))||~aresame(size(x1),size(x3))
    error('The input arrays X1, X2, X3 must all be the same size.')
end

if ~aresame(size(t1),size(t2))||~aresame(size(t1),size(t3))
    error('The input arrays T1, T2, T3 must all be the same size.')
end

maxlent=max([length(t1(:)) length(t2(:)) length(t3(:))]);
maxlenx=max([length(x1(:)) length(x2(:)) length(x3(:))]);

if nargin==7
    if maxlent>1||maxlenx>1
        if ~aresame(size(t1),size(t))&&~aresame(size(x1),size(t))&&~isscalar(t)
            error('The time T must either be the same size as the other input arguments, or a scalar.')
        end
    end
end


numa=x1.*(t2-t3)+x2.*(t3-t1)+x3.*(t1-t2);
denom=(t1-t2).*(t1-t3).*(t2-t3);
a=frac(numa,denom);

numb=x1.*(t2.^2-t3.^2)+x2.*(t3.^2-t1.^2)+x3.*(t1.^2-t2.^2);
b=frac(-numb,denom);

%numc=x1.*(t2.^2.*t3-t3.^2.*t2)+x2.*(t3.^2.*t1-t1.^2.*t3)+x3.*(t1.^2.*t2-t2.^2.*t1);
numc=x1.*t2.*t3.*(t2-t3)+x2.*t3.*t1.*(t3-t1)+x3.*t1.*t2.*(t1-t2);
c=frac(numc,denom);


if nargin==6
    t=-frac(b,2*a);     
end                                   
        
%c=x1-a.*t1.^2-b.*t1;

x=a.*t.^2+b.*t+c;


function[]=quadinterp_test

quadinterp_test_real;
quadinterp_test_complex;
quadinterp_test_extreme;

function[]=quadinterp_test_real


a=1;
b=7;
c=-23;

t1=2;
t2=4;
t3=8;


%t=(-20:.1:20)';
t=[t1 t2 t3];
x=a.*t.^2+b.*t+c;

N=100;
bool=false(N,1);
for i=1:100
   ti=randn(1)*10;
   xi=quadinterp(t(1),t(2),t(3),x(1),x(2),x(3),ti);
   xi2=a.*ti.^2+b.*ti+c;
   bool(i)=aresame(xi,xi2,1e-10);
end

reporttest('QUADINTERP real scalar case',all(bool(i)))


function[]=quadinterp_test_complex
a=1+sqrt(-1)*4;
b=7+sqrt(-1)*2;
c=-23+sqrt(-1)*12;

t1=2;
t2=4;
t3=8;

%t=(-20:.1:20)';
t=[t1 t2 t3];
x=a.*t.^2+b.*t+c;

N=100;
bool=false(N,1);
for i=1:100
   ti=randn(1)*10;
   xi=quadinterp(t(1),t(2),t(3),x(1),x(2),x(3),ti);
   xi2=a.*ti.^2+b.*ti+c;
   bool(i)=aresame(xi,xi2,1e-10);
end

reporttest('QUADINTERP complex scalar case',all(bool(i)))


function[]=quadinterp_test_extreme


a=1;
b=7;
c=-23;

t1=2;
t2=4;
t3=8;


t=[t1 t2 t3];
x=a.*t.^2+b.*t+c;

[xe,te]=quadinterp(t(1),t(2),t(3),x(1),x(2),x(3));

t=(-20:.01:20)';
x=a.*t.^2+b.*t+c;
[xm,im]=min(x);  %This happens to be a minimum  

reporttest('QUADINTERP extrema-finding case',aresame(t(im),te,0.01)&aresame(xm,xe,.01))






