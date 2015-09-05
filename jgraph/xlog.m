function[]=xlog
% XLOG   Sets x-axis scale to logarithmic.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details  

set(gca,'xscale','log')

%Update for compatibility with ZOOL
if isappdata(gcf,'xlimitslock')
   if getappdata(gcf,'xlimitslock')
      h=gca;
      hall=axeshandles;
      for i=1:length(hall)
         axes(hall(i))
         set(gca,'xscale','log')
      end
      axes(h)
   end
end
