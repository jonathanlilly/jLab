function[varargout] = vsum(varargin)
%VSUM  Sum over non-NaN elements along a specified dimension.
%
%   Y=VSUM(X,DIM) takes the sum of all non-NaN elements of X along        
%   dimension DIM. 
%                                                                         
%   [Y,NUM]=VSUM(X,DIM) also outputs the number of non-NaN data points NUM,  
%   which has the same dimension as Y.                             
%
%   [Y1,Y2,...YN]=VSUM(X1,X2,...XN,DIM) also works.
%
%   VSUM(X1,X2,...XN,DIM);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2015 J.M. Lilly --- type 'help jlab_license' for details  
  
if strcmpi(varargin{1}, '--t')
  vsum_test,return
end

dim=varargin{end};

for i=1:length(varargin)-1
  [varargout{i},numi{i}]=vsum1(varargin{i},dim);
end

for i=length(varargin):nargout
  varargout{i}=numi{i-length(varargin)+1};
end

eval(to_overwrite(nargin-1))

function[y,num]=vsum1(x,dim)
x=vswap(x,inf,nan);
y=sum(x,dim,'omitnan');
if nargout==2
    num=sum(isfinite(x),dim,'omitnan');
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=vsum_test
x1=[1 2 ; nan 4];
x2=[nan 6; nan 5];
ans1=[3 4]';
ans2=[6 5]';

vsum(x1,x2,2);
reporttest('VSUM output overwrite', aresame(x1,ans1) && aresame(x2,ans2))

x1=[1 2 ; nan 4];
ans1=[3 4]';
ans2=[2 1]';

[y1,y2]=vsum(x1,2);
reporttest('VSUM sum & num', aresame(y1,ans1) && aresame(y2,ans2))


x1=[1 2 ; 0 4];
ans1=[3 4]';
ans2=[2 2]';

[y1,y2]=vsum(x1,2);
reporttest('VSUM sum & num, no NaNs', aresame(y1,ans1) && aresame(y2,ans2))








