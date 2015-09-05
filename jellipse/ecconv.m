function[varargout]=ecconv(varargin)
%ECCONV  Converts between eccentricity measures.
%  
%   Y=ECCONV(X,STR) converts the eccentricity measure X into eccentricity
%   measure Y as specified by the string STR. 
%  
%   STR is of the form 'in2out' where the names 'in' and 'out' may be
%   chosen from the following table.
%
%       STR        Name                    Definition  
%       ----------------------------------------------------------------
%      'ecc'       Eccentricity            sgn(b)*sqrt(1-b^2/a^2)
%      'ell'       Ellipticity             b/a
%      'lin'       Linearity               sgn(b)*(a^2-b^2)/(a^2+b^2)
%      'circ'      Circularity             2ab/(a^2+b^2) 
%      'nu'        Eccentricity angle      asin(b/sqrt(a^2 + b^2))
%      'rot'       Rotary ratio            sgn(b)*(a-abs(b))/(a+abs(b))
% 
%   where "a" and "b" are the semi-major and semi-minor axes respective. 
%   For example, ECC=ECCONV(L,'lin2ecc') converts the ellipse linearity L
%   into eccentricity ECC.  
%
%   All quantities are positive for mathematically positive rotation
%   and negative for mathematically negative rotation.  Note the linearity
%   is undefined for the case of purely circular rotation a=abs(b);
%
%   For further details, see Lilly and Gascard (2006).
%   ____________________________________________________________________
%   
%   Cell array input/output
%
%   If ECCONV is given cell array input, it returns cell array output.
%
%   Thus X may be a cell arrays of numerical arrays, and the output Y will
%   also be a cell array having an identical size as X.
%   ____________________________________________________________________
%
%   See also ELLCONV, ELLDIFF, TRANSCONV. 
%
%   Usage:   lin=ecconv(b./a,'ell2lin');
%            ecc=ecconv(lin,'lin2ecc');
%
%   'ecconv --t' runs a test.
%   'ecconv --f' generates a sample figure.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2015 J.M. Lilly --- type 'help jlab_license' for details        
 
    
if strcmpi(varargin{1}, '--t')
  ecconv_test,return
end
if strcmpi(varargin{1}, '--f')
  type makefigs_ecconv
  makefigs_ecconv;
  return
end
  
%if nargin==2

x=varargin{1};
str=varargin{2};
ecconv_checkstr(str);

if ~iscell(x)
   y=ecconv_one(x,str);
   y(isinf(x))=inf;
else
    for i=1:length(x)
         y{i,1}=ecconv_one(x{i},str);
         y{i}(isinf(x{i}))=inf;
    end
end

varargout{1}=y;


function[y]=ecconv_one(x,str)

sizex=size(x);
x=x(:);


str1=str(1:strfind(str,'2')-1);
str2=str(strfind(str,'2')+1:end);
if aresame(str1,'circ')&&~aresame(str2,'circ')
    lin=sign(x).*sqrt(1-x.^2);
    eval(['y=lin2' str2 '(x);'])
elseif ~aresame(str1,'circ')&&aresame(str2,'circ')
    eval(['lin=' str1 '2lin(x);'])
    y=sign(lin).*sqrt(1-lin.^2);
else
    eval(['y=' str '(x);'])
end
y=reshape(y,sizex);
  
  
% elseif nargin==3
%   x1=varargin{1};
%   x2=varargin{2};
%   str=varargin{3};
%   ecconv_checkstr(str);
%   eval(['[y1,y2]=' str '(x1,x2);'])
%   varargout{1}=y1;
%   varargout{2}=y2;
% end

function[]=ecconv_checkstr(str)
ii=strfind(str,'2');
if isempty(ii)
  error('STR should be of the form ''ecc2lin''.')
end
str1=str(1:ii-1);
str2=str(ii+1:end);

if ~aresame(str1,'nu') && ~aresame(str1,'ecc') && ~aresame(str1,'lin') && ~ ...
      aresame(str1,'rot') && ~aresame(str1,'ell') && ~aresame(str1,'circ')
  error(['Input parameter name ' str1 ' is not supported.'])
end
if ~aresame(str2,'nu') && ~aresame(str2,'ecc') && ~aresame(str2,'lin') && ~ ...
      aresame(str2,'rot') && ~aresame(str2,'ell') &&~aresame(str2,'circ')
  error(['Output parameter name ' str2 ' is not supported.'])
end

function[nu]=ell2nu(ell)
nu=atan(ell);
function[ell]=nu2ell(nu)
ell=tan(nu);
function[ecc]=ell2ecc(ell)
ecc=nu2ecc(ell2nu(ell));
function[ell]=ecc2ell(ecc)
ell=nu2ell(ecc2nu(ecc));
function[rot]=ell2rot(ell)
rot=nu2rot(ell2nu(ell));
function[ell]=rot2ell(rot)
ell=nu2ell(rot2nu(rot));
function[lin]=ell2lin(ell)
lin=nu2lin(ell2nu(ell));
function[ell]=lin2ell(lin)
ell=nu2ell(lin2nu(lin));

function[ecc]=nu2ecc(nu)
ecc=sign(nu).*sqrt(1-squared(tan(nu)));
vindexinto(ecc,1,find(nu==0),1);
vindexinto(ecc,0,find(abs(nu)==pi/4),1);

function[nu]=ecc2nu(ecc)
nu=sign(ecc).*atan(sqrt(1-squared(ecc)));
vindexinto(nu,pi/4,find(ecc==0),1);

function[lambda]=nu2lin(nu)
lambda=sign(nu).*cos(2*nu);
vindexinto(lambda,1,find(nu==0),1);
vindexinto(lambda,0,find(abs(nu)==pi/4),1);

function[nu]=lin2nu(lambda)
nu=sign(lambda).*atan(sqrt(frac(1-abs(lambda),1+abs(lambda))));
vindexinto(nu,pi/4,find(lambda==0),1);

function[alpha]=nu2rot(nu)
alpha=sign(nu).*frac(1-abs(tan(nu)),1+abs(tan(nu)));
vindexinto(alpha,1,find(nu==0),1);
vindexinto(alpha,0,find(abs(nu)==pi/4),1);

function[nu]=rot2nu(alpha)
nu=sign(alpha).*atan(frac(1-abs(alpha),1+abs(alpha)));
vindexinto(nu,pi/4,find(alpha==0),1);

function[ecc]=lin2ecc(lambda)
ecc=sign(lambda).*sqrt(frac(2*abs(lambda),1+abs(lambda)));
vindexinto(ecc,0,find(lambda==0),1);

function[lambda]=ecc2lin(ecc)
lambda=sign(ecc).*frac(squared(ecc),2-squared(ecc));
vindexinto(lambda,0,find(ecc==0),1);

function[alpha]=lin2rot(lambda)
alpha=nu2rot(lin2nu(lambda));
vindexinto(alpha,0,find(lambda==0),1);

function[lambda]=rot2lin(alpha)
lambda=nu2lin(rot2nu(alpha));
vindexinto(lambda,0,find(alpha==0),1);

function[alpha]=ecc2rot(ecc)
alpha=nu2rot(ecc2nu(ecc));
vindexinto(alpha,0,find(ecc==0),1);

function[ecc]=rot2ecc(alpha)
ecc=nu2ecc(rot2nu(alpha));
vindexinto(ecc,0,find(alpha==0),1);
  
function[x]=rot2rot(x)
function[x]=nu2nu(x)
function[x]=ecc2ecc(x)
function[x]=lin2lin(x)
function[x]=ell2ell(x)
function[x]=circ2circ(x)

 
%/********************************************************
function[]=ecconv_test
str{1}='nu';
str{2}='ecc';
str{3}='lin';
str{4}='rot';
str{5}='ell';

x=rand(100,1)*2-1;
x(end+1,1)=1;  %Add an exact one
x(end+1,1)=0;  %Add an exact zero

clear y z x2 bool

%Turn lambda into all others
for i=1:length(str)
  y(:,i)=ecconv(x,['lin2' str{i}]);
end

%Turn all others into everything else
for i=1:length(str)
  for j=1:length(str)
     z(:,i,j)=ecconv(y(:,i),[str{i} '2' str{j}]);
  end
end

%Turn everything back into lambda
for i=1:length(str)
  for j=1:length(str)
     x2(:,i,j)=ecconv(z(:,i,j),[str{j} '2lin']);
  end
end

for i=1:size(x2,2)
  for j=1:size(x2,2)
    bool(i,j)=aresame(squeeze(x2(:,i,j)),x,1e-14);
  end
end

reporttest('ECCONV', allall(bool))
%\********************************************************


% %Not currently used
% figure
% lambda=(0:.001:1)';
% 
% ecc=ecconv(lambda,'lin2ecc');
% alpha=ecconv(lambda,'lin2rot');
% ecc1=sqrt(2)*sqrt(lambda);
% alpha1=lambda/2;
% 
% plot(lambda,[ecc,alpha,ecc1,alpha1]);
% linestyle k 2k k-- 2k--
% axis([0 1 0 1]),axis square 
% text(0.4,0.7,'\epsilon')
% text(0.7,0.5,'\alpha')
% title('Eccentricity measures')
% xlabel('Ellipse parameter \lambda')
% ylabel('Eccentricity \epsilon or rotary ratio \alpha')
% set(gcf,'paperposition', [2 2 3.5 3.5])
% xtick(.1),ytick(.1),fixlabels(-1)
% fontsize 14 14 14 14

% function[kappa,nu]=ab2kn(a,b)
% kappa=frac(1,sqrt(2))*sqrt(squared(a)+squared(b));
% nu=atan(frac(b,a));
% 
% function[a,b]=kn2ab(kappa,nu)
% a=cos(nu).*kappa*sqrt(2);
% b=sin(nu).*kappa*sqrt(2);
% 
% function[nu]=ab2nu(a,b)
% nu=atan(frac(b,a));

