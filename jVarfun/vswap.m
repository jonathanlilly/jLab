function[varargout]=vswap(varargin)
%VSWAP  VSWAP(X,A,B) replaces A with B in numeric array X.
%
%   VSWAP(X,A,B) replaces A with B in numeric array X.  A and B may be
%   numbers, NAN, +/- INF, NAN+SQRT(-1)*NAN, or INF+SQRT(-1)*INF.
%
%   Note that X may also be a cell array of numeric arrays. 
%
%   [Y1,Y2,...YN]=VSWAP(X1,X2,...XN,A,B) also works.
%
%   VSWAP(X1,X2,...XN,A,B); with no output arguments overwrites the 
%   original input variables.    
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2015 J.M. Lilly --- type 'help jlab_license' for details  

if strcmpi(varargin{1}, '--t')
  vswap_test,return
end
 
  
a=varargin{end-1};
b=varargin{end};

for i=1:length(varargin)-2
  x=varargin{i};
  varargout{i}=swapnum_loop(x,a,b);
end

eval(to_overwrite(nargin-2))  

function[x]=swapnum_loop(x,a,b)
if iscell(x)
    for i=1:length(x)
        x{i}=swapnum1(x{i},a,b);
    end
else
    x=swapnum1(x,a,b);
end

function[x]=swapnum1(x,a,b)
    
    
if isfinite(a)
    bool=(x==a);
else
  if isnan(a)
       if isnan(imag(a))
            bool=isnan(real(x))&isnan(imag(x));
       else
            bool=isnan(x);
       end
  elseif isinf(a)
       if isinf(imag(a))
          bool=isinf(real(x))&isinf(imag(x));
       else
          if a>0
              bool=(isinf(x)&x>0);
          else
              bool=(isinf(x)&x<0);
          end
       end
   end
end


%Convert booleans to double
x=double(x);
x(bool)=b;

function[]=vswap_test
x=(1:10);
ans1=[2 (2:10)];
reporttest('VSWAP num case', aresame(vswap(x,1,2),ans1))

x=[nan (1:10)];
ans1=(0:10);
reporttest('VSWAP nan case', aresame(vswap(x,nan,0),ans1))

x=[nan*(1+sqrt(-1)) nan inf 1:10];
ans1=[0 nan inf 1:10];
reporttest('VSWAP complex nan case', aresame(vswap(x,nan+sqrt(-1)*nan,0),ans1))

x=[inf*(1+sqrt(-1)) nan inf 1:10];
ans1=[0 nan inf 1:10];
reporttest('VSWAP complex inf case', aresame(vswap(x,inf+sqrt(-1)*inf,0),ans1))

x=[nan*(1+sqrt(-1)) nan -inf 1:10];
ans1=[nan*(1+sqrt(-1)) nan 0 1:10];
reporttest('VSWAP negative inf case', aresame(vswap(x,-inf,0),ans1))

