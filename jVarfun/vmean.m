function[varargout] = vmean(varargin)
%VMEAN  Mean over non-NaN elements along a specified dimension.
%
%   Y=VMEAN(X,DIM) takes the mean of all non-NaN elements of X along      
%   dimension DIM. 
%                                                                         
%   [Y,NUM]=VMEAN(X,DIM) also outputs the number of non-NaN data points 
%   NUM, which has the same dimension as X.              
%
%   [Y1,Y2,...YN]=VMEAN(X1,X2,...XN,DIM) or also works, where all the XN
%   are the same size.  
%
%   VMEAN(X1,X2,...XN,DIM);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________
%
%   Weighted means
%
%   VMEAN can also form a weighted mean.
%
%   Y=VMEAN(X,DIM,W) forms the weighted mean of X using weights W, an array 
%   of the same size as X. 
%
%   In this case [Y,NUM]=VMEAN(X,DIM,W) returns total weight used in 
%   forming the mean at each point, rather than the number of data points.
%  
%   Thus VMEAN(X,DIM,W) is the same as VMEAN(X.*W,DIM)./VMEAN(W,DIM).
%
%   [Y1,Y2,...YN]=VMEAN(X1,X2,...XN,DIM,W) also works. 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2015 J.M. Lilly --- type 'help jlab_license' for details    

if strcmpi(varargin{1}, '--t')
  vmean_test;powermean_test;return
end

if ~isscalar(varargin{end})
    w=varargin{end};
    varargin=varargin(1:end-1);
else
    w=[];
end

dim=varargin{end};

for i=1:length(varargin)-1
  if isempty(w)
      [varargout{i},numi{i}]=vsum(varargin{i},dim);
      numi{i}(numi{i}==0)=nan;
      varargout{i}=varargout{i}./numi{i};
  else
      varargout{i}=vsum(varargin{i}.*w,dim);
      numi{i}=vsum(w,dim);
      varargout{i}=varargout{i}./numi{i};
  end
end

for i=length(varargin):nargout
  varargout{i}=numi{i-length(varargin)+1};
end

eval(to_overwrite(nargin-1))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=vmean_test
x1=[1 2 ; nan 4];
x2=[inf 6; nan 5];
ans1=[3/2 4]';
ans2=[6 5]';

vmean(x1,x2,2);
reporttest('VMEAN output overwrite', aresame(x1,ans1) && aresame(x2,ans2))

x1=[1 2 ; nan 4];
ans1=[3/2 4]';
ans2=[2 1]';

[y1,y2]=vmean(x1,2);
reporttest('VMEAN mean & num', aresame(y1,ans1) && aresame(y2,ans2))


function[varargout]=powermean(varargin)
%POWERMEAN  Power-weighted mean along a specified dimension.
%
%   FM=POWERMEAN(F,X,DIM) takes the mean of all finite elements of F along
%   dimension DIM, weighted by the squared magnitude of X.
%   
%   The power-weighted mean of F is defined as 
%
%        FM = SUM (ABS(X)^2.*F, DIM) / SUM(ABS(X)^2, DIM)
%
%   where X and F are arrays of the same size.  
%                                                                         
%   [FM1,FM2,...FMN]=POWERMEAN(F1,F2,...FN,X,DIM) also works.
%
%   POWERMEAN(F1,F2,...FN,X,DIM);  with no output arguments overwrites the 
%   original input variables.
%
%   POWERMEAN with X a set of analytic signals is used to construct the
%   joint instantaneous moments.  For details on these quantities, see
%
%       Lilly and Olhede (2010),  Bivariate instantaneous frequency and
%           bandwidth.  IEEE Trans. Sig. Proc. 58 (2), 591--603.
%
%   See also VMEAN.
%
%   Usage: fm=powermean(f,x,dim);
%          [fm1,fm2]=powermean(f1,f2,x,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    powermean_test,return
end

dim=varargin{end};
x=varargin{end-1};

power=vsum(abs(x).^2,dim);
vswap(power,0,nan);
for i=1:length(varargin)-2
  varargout{i}=vsum(abs(x).^2.*varargin{i},dim)./power;
end

eval(to_overwrite(nargin-2))
 
function[]=powermean_test

bool(1)=aresame(powermean([1 2],[2 3],2),frac(4+2*9,13));
bool(2)=aresame(powermean([1 2],[2 3],2),vmean([1 2],2,[2 3].^2));

reporttest('VMEAN weighted mean',allall(bool))







