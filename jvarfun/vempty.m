function[varargout]=vempty(varargin)
%VEMPTY   Initializes multiple variables to empty sets or empty cell arrays.
%
%   [X1,X2, ... XN]=VEMPTY  is equivalent to
%		
%      X1=[]; X2=[]; .... XN=[];
%
%   thus initializing all the output variables to empty sets.
%
%   [X1,X2, ... XN]=VEMPTY('cell') instead initializes the output variables
%   to empty cell arrays.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004--2013 J.M. Lilly --- type 'help jlab_license' for details

  
if nargin~=0
  if strcmpi(varargin{1}, '--t')
   vempty_test,return
  end
end

str='nocell';
if nargin~=0
    if ischar(varargin{end})
        str=varargin{end};
    end
end

for i=1:nargout
  if strcmpi(str(1:3),'cel')
        varargout{i}=cell(0,1);
  else
        varargout{i}=[];
  end
end

function[]=vempty_test

z=[];
[x,y]=vempty;
reporttest('VEMPTY', all(x==z&y==z))
  
