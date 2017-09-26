function[]=flipy
% FLIPY    Flips the direction of the y-axis
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  

if strcmpi(get(gca,'ydir'),'normal')
  set(gca,'ydir','reverse')
  breverse=1;
elseif strcmpi(get(gca,'ydir'),'reverse')
  set(gca,'ydir','normal')
  breverse=0;
end

%Update for compatibility with ZOOL
if isappdata(gcf,'ylimitslock')
   if getappdata(gcf,'ylimitslock')
      h=gca;
      hall=axeshandles;
      for i=1:length(hall)
         if breverse	 
             set(hall(i),'ydir','reverse')
         else
              set(hall(i),'ydir','normal')
         end
      end
   end
end
