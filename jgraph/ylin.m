function[]=ylin
% YLIN   Sets y-axis scale to linear.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details  
    
  
set(gca,'yscale','linear')

%Update for compatibility with ZOOL
if isappdata(gcf,'ylimitslock')
   if getappdata(gcf,'ylimitslock')
      h=gca;
      hall=axeshandles;
      for i=1:length(hall)
         axes(hall(i))
         set(gca,'yscale','linear')
      end
      axes(h)
   end
end
