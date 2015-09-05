function[varargout]=vcolon(varargin)
%VCOLON  Condenses its arguments, like X(:).
%
%   [Y1,Y2, ... YN]=VCOLON(X1,X2, ... XN) is equivalent to   
%		
%      Y1=X1(:); Y2=X2(:); ... YN=XN(:);    
%
%   VCOLON(X1,X2,...XN) with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details    
    
if strcmpi(varargin{1}, '--t')
   vcolon_test,return
end
 
for i=1:nargin
  x=varargin{i};
  varargout{i}=x(:);
end

eval(to_overwrite(nargin));


function[]=vcolon_test
x=[1 3; 2 4];
y=x;

vcolon(x,y);
reporttest('VCOLON ', all(x==(1:4)'&y==(1:4)'))
