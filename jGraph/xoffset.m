function[]=xoffset(i1,i2)
%XOFFSET  Offsets lines in the x-direction after plotting.
%	    	    	    
%   XOFFSET allows data to be manipulated after it has been plotted, by
%   operating on the data stored in the figure itself.  This may be used 
%   only with 2-D line plots (not symbol plots, contour plots, etc.). 
%
%   XOFFSET offsets all lines in the current axes a specified amount.
%		
%   XOFFSET(N) or XOFFSET N, where N is a *real* number, offsets each
%   column of the X-data by an amount N from the previous column.  N must 
%   be a number, not a variable whose value is a number.
%
%   XOFFSET(N) or XOFFSET N, where N is an *imaginary* number, offsets
%   each column of the X-data by an amount imag(N)*DX from the previous 
%   column, where DX is the X-axis length of the original (unoffset) plot.
%
%   XOFFSET(M,N) or XOFFSET M N, applies the same offsets to each of M 
%   groups of lines. For example, for complex-valued data Z use UVPLOT(Z) 
%   followed by XOFFSET 2 N.
%
%   XOFFSET with no arguments returns the data to its original (unoffset) 
%   orientation [as do XOFFSET(0) and XOFFSET(0i)].
%   
%   XOFFSET also allow multiple axes to be manipulated simultaneously.
%
%   XOFFSET LOCK and XOFFSET UNLOCK lock and unlock all Y-data in the
%   current figure. When LOCK is on, calls to XOFFSET are applied to each
%   subplot individually.
%
%   See also YOFFSET, LINESTYLE.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details  

if strcmpi(i1,'--t')
    return
end

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

bcontinue=1;
if nargin>0
   if ischar(i1)
      if isempty(str2double(i1))
           i1=deblank(i1);
           i1=fliplr(deblank(fliplr(i1)));	
           if strcmpi(i1,'lock')
                if isappdata(gca,'lastxoffset')
                    lastoff=getappdata(gca,'lastxoffset');
                    setappdata(gcf,'xoffsetlock',1)
                    xoffsetloop(0,N);
                    xoffsetloop(lastoff,N); 
                else
                    setappdata(gcf,'xoffsetlock',1)
                end		 
                bcontinue=0;
            end
            if strcmpi(i1,'unlock')
                setappdata(gcf,'xoffsetlock',0)
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
   %If the plot has been changed, the old "lastxoffset" may no
   %longer be correct.  If so, we set it to zero
       setappdata(gca,'lastxoffset',0)		 
   end
   xoffsetloop(delta,N);
end




function[]=xoffsetloop(delta,N)
h=gca;
if isappdata(gcf,'xoffsetlock')
   if getappdata(gcf,'xoffsetlock')
      h=axeshandles(gcf);	
   end 
end

for i=1:length(h)
    axes(h(i))
    if isappdata(gca,'lastxoffset')
       lastdelta=getappdata(gca,'lastxoffset');
       xoffsetapply(-1*lastdelta,N);
       setappdata(gca,'lastxoffset',0)
    end
    xoffsetapply(delta,N);
end

axes(h(1))

function[]=xoffsetapply(delta,N)
h=linehandles;

if ischar(delta),delta=str2double(delta);end
xx=get(h,{'xdata'});

if length(delta)==1 && ~isreal(delta)
   ax=axis;
   dx=ax(2)-ax(1);
   delta=dx.*imag(delta)/100;
end

if length(delta)==1
   for i=length(xx):-1:1
       xx{i}=xx{i}+delta*(mod(length(xx)-i,round(length(xx)/N)));
   end
elseif length(delta)==length(h)
   for i=length(xx):-1:1
       xx{i}=xx{i}+delta(i);
   end
else
   error('Number of line handles does not equal offset array length.')
end

set(h,{'xdata'},xx);
setappdata(gca,'lastxoffset',delta)		
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
  


