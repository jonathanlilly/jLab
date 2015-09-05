function[h]=hlines(y,sty)
%HLINES   Add horizontal lines to a plot.
%
%   HLINES(Y,STYLE) puts horizontal lines at locations Y having style
%   described by the string STYLE.  HLINES draws the lines to fit the
%   current axis.
%
%   The SYTLE strings follow the format specified in LINESTYLE.  Thus
%   STYLE='2b--' draws blue dotted lines of width 2.
%
%   HLINES(Y) uses the default STYLE='g--', a green dashed line.
%
%   H=HLINES(...) returns an array of handles H to the lines.
%
%   See also VLINES, DLINES.	
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2008 J.M. Lilly --- type 'help jlab_license' for details  

if strcmpi(y,'--t')
    return
end

if nargin==1
	sty='g--';
end

hold on
y=y(:)';
ax=axis;
x1=ax(1)*ones(size(y));
x2=ax(2)*ones(size(y));
y=[y;y];
x=[x1;x2];
h=plot(x,y);
linestyle(h,sty);
axis(ax)

if nargout==0
   clear h
end

