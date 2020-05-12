function[varargout] = vstd(varargin)
%VSTD Standard deviation over non-NaN elements along a specfied dimension.
%
%   Y=VSTD(X,DIM) takes the standard deviation of all non-NaN elements of X 
%   along dimension DIM.
%                                                                         
%   [Y,NUM]=VSTD(X,DIM) also outputs the number of non-NaN data points NUM, 
%   which is the same size as X.
%                                                                         
%   VSTD uses the "1/N" normalization, where N is the number of data points,
%   rather than the "1/(N-1)" normalization.
%
%   If X is complex-valued, X=XR+1i*XI, the stardard deviation as it is
%   formally defined is not particularly meaningful.  Therefore, the output
%   of VSTD is redefined such that Y=VSTD(XR,DIM)+1i*VSTD(XI,DIM).
%
%   [Y1,Y2,...YN]=VSTD(X1,X2,...XN,DIM) also works.
%
%   VSTD(X1,X2,...XN,DIM); with no arguments overwrites the original 
%   input variables.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2020 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmpi(varargin{1}, '--t')
  vstd_test,return
end

dim=varargin{end};

for i=1:length(varargin)-1
   if isreal(varargin{i})
       [varargout{i},numi{i}]=vmoment(varargin{i},2,dim);
       varargout{i}=sqrt(varargout{i});
   else
       [sigx2,~]=vmoment(real(varargin{i}),2,dim);
       [sigy2,numi{i}]=vmoment(imag(varargin{i}),2,dim);
       varargout{i}=sqrt(sigx2)+1i*sqrt(sigy2);
   end
   if isreal(varargin{i})
       varargout{i}=real(varargout{i});
   end
end

for i=length(varargin):nargout
  varargout{i}=numi{i-length(varargin)+1};
end

eval(to_overwrite(nargin-1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=vstd_test
x1=[1 2 3 nan];
x2=x1;
ans1=sqrt(2/3);
vstd(x1,x2,2);
reporttest('VSTD output overwrite', aresame(x1,ans1) && aresame(x2,ans1))
ans2=3;

x1=[1 2 3 nan];
[y1,y2]=vstd(x1,2);
reporttest('VSTD std & num', aresame(y1,ans1) && aresame(y2,ans2))

cv=randn(100,4)+sqrt(-1)*randn(100,4);
clear cvstd
cvstd(:,1)=conj(vstd(cv,1))';
cvstd(:,2)=vstd(real(cv),1)';
cvstd(:,3)=vstd(imag(cv),1)';
%cvstd(:,4)=sqrt(cvstd(:,2).^2+cvstd(:,3).^2);
cvstd(:,4)=cvstd(:,2)+1i*cvstd(:,3);
bool=aresame(cvstd(:,1),cvstd(:,4),1e-10);
reporttest('VSTD complex-valued data', bool)
