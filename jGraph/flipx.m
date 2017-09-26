function[]=flipx
% FLIPX    Flips the direction of the x-axis
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  

if strcmpi(get(gca,'xdir'),'normal')
  set(gca,'xdir','reverse')
  breverse=1;
elseif strcmpi(get(gca,'xdir'),'reverse')
  set(gca,'xdir','normal')
  breverse=0;
end

%Update for compatibility with ZOOL
if isappdata(gcf,'xlimitslock')
   if getappdata(gcf,'xlimitslock')
      h=gca;
      hall=axeshandles;
      for i=1:length(hall)
         if breverse	 
             set(hall(i),'xdir','reverse')
         else
             set(hall(i),'xdir','normal')
         end
      end
   end
end
