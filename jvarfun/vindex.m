function[varargout]=vindex(varargin)
%VINDEX  Indexes an N-D array along a specified dimension.
%
%   Y=VINDEX(X,INDEX,DIM) indexes the multidimensional array X along     
%   dimension DIM. This is equivalent to   
%		
%		    1 2       DIM     DIMS(X)
%		    | |        |         |
%		Y=X(:,:, ... INDEX, ..., :);		
%
%   where the location of INDEX is specified by DIM.
%
%   VINDEX is defined to return an empty array if INDEX is empty.
%  
%   Note that VINDEX does not index along singleton dimensions, thus
%   when X is a column vector, VINDEX(X,INDEX,2) returns X.  
%
%   [Y1,Y2,...YN]=VINDEX(X1,X2,...XN,INDEX,DIM) also works.
%
%   VINDEX(X1,X2,...XN,INDEX,DIM); with no output arguments overwrites 
%   the original input variables.
%
%   VINDEX also supports logical indexing with INDEX a boolean array
%   of the same size as the dimension being indexed.
%
%   See also VINDEXINTO, SQUEEZE, DIMS, PERMUTE, SHIFTDIM.
%
%   'vindex --t' runs a test.
%
%   Usage:  y=vindex(x,index,dim);
%           [y1,y2,y3]=vindex(x1,x2,x3,index,dim);
%           vindex(x1,x2,x3,index,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2014 J.M. Lilly --- type 'help jlab_license' for details    


if strcmpi(varargin{1}, '--t')
  vindex_test,return
end

index=varargin{end-1};
dim=varargin{end};

nvars=nargin-2;
vars=varargin(1:nvars);

%eval(to_grab_from_caller(2))  %assigns vars, varnames, nvars

for i=1:nvars
  varargout{i}=vindex1(vars{i},index,dim);
end

eval(to_overwrite(nargin-2));

%now to_overwrite uses
%varnames not inputnames ... and does not take input argument
   
function[y]=vindex1(x,index,dim)  

y=[];
if ~isempty(x)
    if size(x,dim)==1
        %X has only one element along dimension DIM; do nothing
        if ~isempty(index)
            y=x;
        end
    else
        %You would think Matlab would provide a simpler way to do this.
        str='y=x(';
        ndx=length(find(size(x)>=1));
        if ~isempty(index)
          for i=1:ndx
              if i~=dim
                  str=[str ':,'];
              else
                  str=[str 'index,'];
              end
          end              
          str=[str(1:end-1) ');'];
          eval(str);
        end
    end
end

function[]=vindex_test
x1=[1 1; 2 2];
index=2;
ans1=x1(:,index);
vindex(x1,index,2);
reporttest('VINDEX col case', aresame(x1,ans1))

x1=[1 1; 2 2];
index=2;
ans1=x1(index,:);
vindex(x1,index,1);
reporttest('VINDEX row case', aresame(x1,ans1))

x1=[1 1; 2 2];
index=[false; true];
ans1=x1(:,index);
vindex(x1,index,2);
reporttest('VINDEX col case, logical', aresame(x1,ans1))

x1=[1 1; 2 2];
index=[false; true];
ans1=x1(index,:);
vindex(x1,index,1);
reporttest('VINDEX row case, logical ', aresame(x1,ans1))

x1=[1 2];
index=1:10;
ans1=x1;
vindex(x1,index,1);
reporttest('VINDEX row vector indexed along rows case', aresame(x1,ans1))

x1=[1 2]';
index=1:10;
ans1=x1;
vindex(x1,index,2);
reporttest('VINDEX column vector indexed along columns case', aresame(x1,ans1))

x1=[1 2]';
index=[];
ans1=[];
vindex(x1,index,2);
reporttest('VINDEX empty index case', aresame(x1,ans1))

x1=[];
index=(1:2);
ans1=[];
vindex(x1,index,2);
reporttest('VINDEX empty array case', aresame(x1,ans1))
