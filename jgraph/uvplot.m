function[h1,h2]=uvplot(varargin)
%UVPLOT  Plots the real and imaginary parts of a signal on the same axis.
%
%   UVPLOT(T,CV) plots the real and imaginary components of the
%   complex-valued data CV=U+iV versus time T.  Here T is a column
%   vector and CV is a matrix with LENGTH(T) rows.
%
%   UVPLOT(T,CV,STY1,STY2) optionally plots the real components of CV,
%   U, using style STY1, and the imaginary components V using STY2.
%   STY1 and STY2 are strings following the format in LINESTYLE.
%
%   UVPLOT(T,CV,STY) uses STY for the linestyles of both lines.
%
%   UVPLOT(T,U,V,...), where U and V are real matrices, also works.  
%
%   [H1,H2]=UVPLOT(...) returns handles to the line plots of U and V.
%
%   As an example,
%
%        load bravo94
%        cv=vfilt(bravo.rcm.cv(:,3),24);
%        num=yf2num(bravo.rcm.yearf)-datenum(1994,1,1);
%        uvplot(num,cv,'2b','r-.'),xlim([120 340])
%
%   returns a color version of Fig.4b of Lilly and Rhines JPO 2001.
% 
%   Usage:  uvplot(t,cv);
%           uvplot(t,cv,sty1,sty2);
%           uvplot(t,cv,sty1,sty2,dy);
%           uvplot(t,u,v,sty1,sty2,dy);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2012 J.M. Lilly --- type 'help jlab_license' for details    

if strcmpi(varargin{1},'--t')
    return;
end
%Depricated
%   UVPLOT(T,CV,STY1,STY2,DY) also offsets the successive columns of
%   CV by an amount DY, with positive DY implying the first column
%   will at the top.  DY may also be a vector specifying an offset for
%   each column.

na=nargin;
boolhold=ishold(gca);

if isrow(varargin{1})
  varargin{1}=varargin{1}(:);
end

offs=0;
if nargin>1
  if ~ischar(varargin{na});
    if length(varargin{na})~=size(varargin{1},1);
       offs=varargin{na};
       na=na-1;
    end
  end
end

%/********************************************************
%Line style arguments
%sty1='k';
%sty2='2G';

%if verLessThan('matlab','8.4.0')
%sty1='b';
%sty2='r';
%sty3='g';

if ischar(varargin{na})
   sty2=varargin{na};
   na=na-1;
   sty1=sty2;
end

if nargin>1
  if ischar(varargin{na})
    sty1=varargin{na};    
    na=na-1;
  end
end
%\********************************************************

if na==1
	u=real(varargin{1});
	v=imag(varargin{1});
	ii=1:length(u);
elseif na==2
	ii=varargin{1};
	u=real(varargin{2});
	v=imag(varargin{2});
elseif na==3
	ii=varargin{1};
	u=varargin{2};
	v=varargin{3};
end

% if size(u,2)==length(offs) || length(offs==1);
%   offs=offs(:);
%   if length(offs)==1
%     offs=flipud(cumsum(offs+zeros(size(u,2),1)));
%   end
%   u=vadd(u,offs',1);
%   v=vadd(v,offs',1);
% end

h1=plot(ii,u); hold on
h2=plot(ii,v);
%linestyle(h1,sty1);
%linestyle(h2,sty2);

if nargout==0
  clear h1 h2
end

	
if ~boolhold
  hold off
end
%h(1)=plot(ii,u,'b');hold on;h(2)=plot(ii,v,'g');
 
