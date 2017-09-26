function[varargout]=vindexinto(varargin)
%VINDEXINTO  Indexes into N-D array along a specified dimension.
%
%   Y=VINDEXINTO(Y,X,INDEX,DIM) indexes the array X into the multi- 
%   dimensional array Y along dimension DIM. This is equivalent to   
%		
%		    1 2       DIM     DIMS(X)
%		    | |        |         |
%		  Y(:,:, ... INDEX, ..., :)=X;		
%
%   where the location of INDEX is specified by DIM.  X must have the
%   exact same size as the block it replaces, or be a scalar.  
%   
%   VINDEXINTO is defined leave Y unchanged if INDEX is empty, and to 
%   ignore NANs and INFs in INDEX.  
%
%   VINDEXINTO(Y,X,INDEX,0) with DIM=0 is equivalent to Y(INDEX)=X;
%  
%   [Y1,Y2,...YN]=VINDEXINTO(Y1,Y2,...YN,X1,X2,...XN,INDEX,DIM) also 
%   works.
%
%   VINDEXINTO(Y1,Y2,...YN,X1,X2,...XN,INDEX,DIM); with no output 
%   arguments overwrites the original input variables ,Y2,...YN.
%
%   See also VINDEX, SQUEEZE, DIMS, PERMUTE, SHIFTDIM
%
%   'vindexinto --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2014 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmpi(varargin{1}, '--t')
  vindexinto_test,return
end

index=varargin{end-1};
dim=varargin{end};

nvars=nargin-2;
if ~iseven(nvars)
  error('There should be an even number of input arguments.')
end

vars=varargin(1:nvars/2);
fromvars=varargin(nvars/2+1:nvars);

for i=1:nvars/2
varargout{i}=vindexinto1(vars{i},fromvars{i},index,dim);
end
eval(to_overwrite(nvars/2));

function[y]=vindexinto1(y,x,index,dim)  
%You would think Matlab would provide a simpler way to do this.

if dim==0
    y(index)=x;
else
    str='y(';

    %[index,sorter]=sort(index);
    %x=x(sorter);

    ii=find(isfinite(index));
    if ~isempty(ii)
      index=index(ii);
    else
      x=[];
      index=[];
    end


    ndx=length(find(size(y)>=1));
    if ~isempty(index) && ~isempty(x)
      for i=1:ndx
          if i~=dim
              str=[str ':,'];
          else
        str=[str 'index,'];
          end
      end
      str=[str(1:end-1) ')=x;'];
      eval(str);
    end
end

% if any(~isfinite(y(:)))
%   if all(isreal(y(isfinite(y))))
%       y(~isfinite(y(:)))=nan;
%   else
%       y(~isfinite(y(:)))=nan+sqrt(-1)*nan;
%   end
% end      
    
% if any(isnan(y(:)))
%   if all(isreal(y(isfinite(y))))
%       y(isnan(y(:)))=nan;
%   else
%       y(isnan(y(:)))=nan+sqrt(-1)*nan;
%   end
%end    

function[]=vindexinto_test
y1=[1 1; 2 2];
index=2;
x=[5 6]';
ans1=y1;
ans1(:,2)=x;
vindexinto(y1,x,index,2);
reporttest('VINDEXINTO col case', aresame(y1,ans1))

y1=[1 1; 2 2];
index=2;
x=[5 6];
ans1=y1;
ans1(2,:)=x;
vindexinto(y1,x,index,1);
reporttest('VINDEXINTO row case', aresame(y1,ans1))

clear y1 ans1
y1(:,:,1)=[1 2; 3 4];
y1(:,:,2)=2*[1 2; 3 4];
index=1;
x=[5 6]';x=vrep(x,2,3);
ans1(:,:,1)=[5 6; 3 4];
ans1(:,:,2)=[5 6; 2*3 2*4];
vindexinto(y1,x,index,1);
reporttest('VINDEXINTO 3-D col case', aresame(y1,ans1))

vindexinto(y1,[],index,1);
reporttest('VINDEXINTO empty x case', aresame(y1,ans1))


y1=[1 1; 2 2];
index=[1 4];
x=[3 5];
ans1=[3 1; 2 5];
vindexinto(y1,x,index,0);
reporttest('VINDEXINTO DIM=0 case', aresame(y1,ans1))
