function[]=yoffset(i1,i2)
%YOFFSET  Offsets lines in the y-direction after plotting.
%	    	    	    
%   YOFFSET allows data to be manipulated after it has been plotted, by
%   operating on the data stored in the figure itself.  This may be used 
%   only with 2-D line plots (not symbol plots, contour plots, etc.). 
%
%   YOFFSET offsets all lines in the current axes a specified amount.
%		
%   YOFFSET(N) or YOFFSET N, where N is a *real* number, offsets each
%   column of the Y-data by an amount N from the previous column.  N must 
%   be a number, not a variable whose value is a number.
%
%   YOFFSET(N) or YOFFSET N, where N is an *imaginary* number, offsets
%   each column of the Y-data by an amount imag(N)*DY from the previous 
%   column, where DY is the Y-axis length of the original (unoffset) plot.
%
%   YOFFSET(M,N) or YOFFSET M N, applies the same offsets to each of M 
%   groups of lines. For example, for complex-valued data Z use UVPLOT(Z) 
%   followed by YOFFSET 2 N.
%
%   YOFFSET with no arguments returns the data to its original (unoffset) 
%   orientation [as do YOFFSET(0) and YOFFSET(0i)].
%   
%   YOFFSET also allow multiple axes to be manipulated simultaneously.
%
%   YOFFSET LOCK and YOFFSET UNLOCK lock and unlock all Y-data in the
%   current figure. When LOCK is on, calls to YOFFSET are applied to each
%   subplot individually.
%
%   See also YOFFSET, LINESTYLE.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details  



if strcmpi(i1,'--t')
    return
end

bcontinue=1;
if nargin==1;
  N=10000;  %Very large number of lines
end
if nargin==2
  N=i1;
  i1=i2;
  if ischar(N)
    N=str2double(N);
  end 
end

if nargin>0
   if ischar(i1)
      if isempty(str2double(i1))
           i1=deblank(i1);
           i1=fliplr(deblank(fliplr(i1)));	
           if strcmpi(i1,'lock')
               if isappdata(gca,'lastyoffset')
                    lastoff=getappdata(gca,'lastyoffset');
                    setappdata(gcf,'yoffsetlock',1)
                    yoffsetloop(0,N);
                    yoffsetloop(lastoff,N); 
               else
                    setappdata(gcf,'yoffsetlock',1)
               end		 
               bcontinue=0;
           end
           if strcmpi(i1,'unlock')
               setappdata(gcf,'yoffsetlock',0)
               bcontinue=0;
           end
      end
   end
end

if bcontinue	
   if nargin>0
      delta=i1;
   else 
      delta=0;
   end

   h=linehandles;
   if plotmodified(h)
   %If the plot has been changed, the old "lastyoffset" may no
   %longer be correct.  If so, we set it to zero
       setappdata(gca,'lastyoffset',0)		 
   end
   yoffsetloop(delta,N);
end



function[]=yoffsetloop(delta,N)

h=gca;
if isappdata(gcf,'yoffsetlock')
   if getappdata(gcf,'yoffsetlock')
      h=axeshandles(gcf);	
   end 
end

for i=1:length(h)
    axes(h(i))
    if isappdata(gca,'lastyoffset')
       lastdelta=getappdata(gca,'lastyoffset');
       yoffsetapply(-1*lastdelta,N);
       setappdata(gca,'lastyoffset',0)
    end
    yoffsetapply(delta,N);
end

axes(h(1))

function[]=yoffsetapply(delta,N)
h=linehandles;
     
if ischar(delta),delta=str2double(delta);end
yy=get(h,{'ydata'});

if length(delta)==1 && ~isreal(delta)
   ax=axis;
   dy=ax(4)-ax(3);
   delta=dy.*imag(delta)/100;
end

if length(delta)==1
   for i=length(yy):-1:1
       yy{i}=yy{i}+delta*(mod(length(yy)-i,round(length(yy)/N)));
   end
elseif length(delta)==length(h)
   for i=length(yy):-1:1
       yy{i}=yy{i}+delta(i);
   end
else
   error('Number of line handles does not equal offset array length.')
end

set(h,{'ydata'},yy);
setappdata(gca,'lastyoffset',delta)		
setappdata(gca,'lasthandles',h)		

function[b]=plotmodified(h)
b=0;  
if isappdata(gca,'lasthandles')
  b=1;
  h2=getappdata(gca,'lasthandles');
  if length(h2)==length(h)
      if all(h2==h)
          b=0;
      end
  end
end
  
