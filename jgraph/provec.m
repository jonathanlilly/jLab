function[data,h,hc]=provec(varargin)
%PROVEC  Generate progressive vector diagrams (simple and fancy).
%
%   Simple provecs:
%
%     PROVEC(U,V,DELTAT) generates a simple progressive vector diagram
%     plotting cumsum(U*DELTAT) vs cumsum(V*DELTAT).  U and V are
%     column vectors, or matrices with time oriented in columns.
%     DELTAT (optional) must be in units of seconds and defaults to
%     one hour (3600 s).  U and V must have units of cm/s.
%  
%     PROVEC(CV,DELTAT), where CV=U+iV, also works.
%
%     INT=PROVEC(...) outputs the integrated dispacement INT in km.
%	 
%   Fancy provecs:
%
%     A fancy provec use SCATTER to plot change the color and/or sizes
%     of the points according to another parameter, say density.  A
%     colorbar is also plotted.
%
%     PROVEC(U,V,C,DELTAT) uses C as the symbol color.  The array C
%     must have the same size as U and V.
%
%     PROVEC(U,V,C,S,DELTAT) also uses S as the symbol size.  The
%     array S must have the same size as U and V.
%	
%     PROVEC(CV,C,DELTAT) and PROVEC(CV,C,S,DELTAT) also work.
%
%     Note that since SCATTER is slow for large datasets, it is useful
%     to decimate the data before plotting.  This is accomplished
%     using PROVEC(...,INDEX,[DELTAT]), where INDEX is placed just
%     before the optional DELTAT.  Then only the points INT(INDEX,:)
%     of the integrated trajectory are plotted.
%
%     [INT,H,HC]=PROVEC returns the intergrated displacement INT, the
%     handle H to the data, and the handle HC to the colorbar axis.
%
%   PROVEC(DCOL,...), where DCOL is a scalar or a row vector,
%   specifies that the columns of the data are to by offset by amount
%   DCOL after plotting, for both simple and fancy provecs.  
%
%   If DCOL is a scalar, then successive columns after the first are
%   offset by amount -DCOL.  If DCOL is a row vector, then it
%   specifies the offset for each column of the data.
%
%   As an example,
%  
%          load bravo94
%          th=100*detrend(bravo.cat.th(:,3));
%          [int,h,hc]=provec(bravo.rcm.cv(:,3),th,20+0*th,[1:10:4000]);
%          caxis([-8 8])
%
%   makes part of Figure 6b of Lilly and Rhines (2002) JPO.
%         
%   Usage: int=provec(cv);
%          [int,h,hc]=provec(cv,c,index);
%          [int,h,hc]=provec(cv,c,s,index);
%          [int,h,hc]=provec(cv,c,s,index,deltat);
%          [int,h,hc]=provec(dcol,cv,c,s,index,deltat);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1999--2012 J.M. Lilly --- type 'help jlab_license' for details        

inunits=100;
outunits=1000;
factor=100*1000;

bcolor=0;
defsize=5;

na=nargin;
deltat=3600;

%/********************************************************
%Look for initial row vector
offs=0;
if isrow(varargin{1}) || isscalar(varargin{1});
  offs=varargin{1};
  na=na-1;
  varargin=varargin(2:end);
end
%\********************************************************

if length(varargin{end})==1
  deltat=varargin{end};
  na=na-1;
end

index=[];
if length(varargin{na})~=length(varargin{1})
  index=varargin{na};
  na=na-1;
end

c=[];
if isreal(varargin{1})
  data=varargin{1}+sqrt(-1)*varargin{2};
  if na>2
      c=varargin{3};
  end
  if na>3
      s=varargin{4};
  else
      s=defsize+zeros(size(data));
  end 
else
  data=varargin{1};
  if na>1
      c=varargin{2};
  end
  if na>2
      s=varargin{3};
  else
      s=defsize+zeros(size(data));
  end 
end

data=vswap(data,nan,0);
data=cumsum(data)*deltat/factor;
if ~(length(offs)==1&&allall(abs(offs)==0))
  if isscalar(offs)
    offs=-(0:1:size(data,2)-1)*offs;
  end
  if length(offs)~=size(data,2)
    error('Length of DCOL must equal number of columns of the data.')
  end
  for i=1:size(data,2)
      data(:,i)=data(:,i)+offs(i);
  end
end


if ~isempty(index);
  data=data(index,:);
  if ~isempty(c);
    c=c(index,:);
  end
  if ~isempty(s);
    s=s(index,:);
  end
end


if isempty(c)
	h=plot(data);
else 
	bcolor=1;
        vcolon(data,s,c);
	if any(isnan(s))
	  error('S cannot contain NANs.')
	end
	if any(isnan(c))
	  error('C cannot contain NANs.')
	end
	h=scatter(real(data),imag(data),s,c,'filled');
end

xlabel('Displacement eastward (km)')
ylabel('Displacement northward (km)')

set(gca,'box','on')
set(gca,'dataaspectratio',[1 1 1])
pos=get(gca,'position');

%put a colorbar if needed
if bcolor
	ax=gca;
	hc=colorbar;
	posc=get(hc,'position');
	set(hc,'position',[posc(1) pos(2) posc(3) pos(4)])
	%the above doesn't take care of the relative size
	%problem--- instead, try to change aspect ratio
	axes(ax)
end

hold on

if nargout ==0
    clear data h hc
end
  
  
  
  
  




