function[h,hc]=hodograph(varargin)
%HODOGRAPH   Generate hodograph plots (simple and fancy).
%
%   A hodograph of a velocity time series is defined as a polar
%   coordinate plot of the eastward component versus the northward.
%  
%   Simple hodographs:
%
%     HODOGRAPH(U,V) generates a simple hodograph plotting U versus V.
%     U and V are column vectors, or matrices with time oriented in
%     columns.  
%  
%     HODOGRAPH(CV), where CV=U+iV, also works.
%
%   Fancy hodographs:
%
%     A fancy hodograph use SCATTER to plot change the color and/or
%     sizes of the points according to another parameter, say density.
%     A colorbar is also plotted.
%
%     HODOGRAPH(U,V,C) uses C as the symbol color.  The array C must
%     have the same size as U and V.
%
%     HODOGRAPH(U,V,C,S) also uses S as the symbol size.  The array S
%     must have the same size as U and V.
%	
%     HODOGRAPH(CV,C) and HODOGRAPH(CV,C,S) also work.
%
%     Note that since SCATTER is slow for large datasets, it is useful
%     to decimate the data before plotting.  This is accomplished
%     using HODOGRAPH(...,INDEX). Then only the points CV(INDEX,:) are
%     plotted.
%
%     [H,HC]=FHODOGRAPH(...) returns the handle H to the data, and the
%     handle HC to the colorbar axis.
%
%   As an example,
%  
%          load bravo94
%          cv=vfilt(bravo.rcm.cv(:,3),24);
%          th=100*detrend(bravo.cat.th(:,3));
%          [h,hc]=hodograph(cv,th,20+0*th,[1:4000]);
%          caxis([-8 8])
%         
%   makes Figure 6a of Lilly and Rhines (2002) JPO.
%         
%   Usage: hodograph(cv);
%          [h,hc]=hodograph(cv,c,index);
%          [h,hc]=hodograph(cv,c,s,index);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2012 J.M. Lilly --- type 'help jlab_license' for details  


  
inunits=100;
outunits=1000;

bcolor=0;
defsize=5;

na=nargin;
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

if ~isempty(index);
  data=data(index,:);
  if ~isempty(c);
    c=c(index,:);
  end
  if ~isempty(s);
    s=s(index,:);
  end
end


u0=max(max(abs(data(:))));
%ceil(u0/5)*5;

bfixlabels=0;
if ~ishold
	h=polar(0,u0);
	hold on
	bfixlabels=1;
else
	h=gca;
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

xlabel('Eastward velocity (km)')
ylabel('Northward velocity (km)')

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
  
  
