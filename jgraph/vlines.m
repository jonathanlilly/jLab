function[h]=vlines(x,sty)
%VLINES   Add vertical lines to a plot.
%
%   VLINES(X,STYLE) puts vertical lines at locations X having style
%   described by the string STYLE.  VLINES draws the lines to fit the
%   current axis.
%
%   The SYTLE strings follow the format specified in LINESTYLE.  Thus
%   STYLE='2b--' draws blue dotted lines of width 2.
%
%   VLINES(X) uses the default STYLE='g--', a green dashed line.
%
%   H=VLINES(...) returns an array of handles H to the lines.
%
%   See also HLINES, DLINES.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2008 J.M. Lilly --- type 'help jlab_license' for details    

if strcmpi(x,'--t')
    return
end

if nargin==1
	sty='g--';
end
hold on

x=x(:)';  
ax=axis;
y1=ax(3)*ones(size(x));
y2=ax(4)*ones(size(x));
x=[x;x];
y=[y1;y2];
h=plot(x,y);
linestyle(h,sty);
axis(ax)

if nargout==0
   clear h
end
