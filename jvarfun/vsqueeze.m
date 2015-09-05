function[varargout]=vsqueeze(varargin)
%VSQUEEZE   Squeezes multiple input arguments simultaneously.
%
%   [Y1,Y2, ... YN]=VSQUEEZE(X1,X2, ... XN) is equivalent to   
%		
%      Y1=SQUEEZE(X1); Y2=SQUEEZE(X2); ... YN=SQUEEZE(XN);
%
%   VSQUEEZE(X1,X2,...XN) with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details    

  
if strcmpi(varargin{1}, '--t')
   vsqueeze_test,return
end
 
for i=1:nargin
  x=varargin{i};
  varargout{i}=squeeze(x);
end

eval(to_overwrite(nargin));


function[]=vsqueeze_test
x=zeros(10,10,10);
y=x;

z=squeeze(x(:,4,:));

vindex(x,y,4,2);
vsqueeze(x,y);
reporttest('VSQUEEZE', all(x==z&y==z))
  
