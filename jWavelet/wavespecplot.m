function[h,hl]=wavespecplot(varargin)
%WAVESPECPLOT  Plot of wavelet spectra together with time series.
%
%   WAVESPECPLOT(T,X,P,W) where T is the time, X is a time series, and W is 
%   the wavelet transform at periods P, makes a two-component plot.  The 
%   upper subplot has the time series X plotted against its time axis T.
%   The lower subplot has the transform W, or its modulus ABS(W) if W is 
%   complex-valued, plotted versus the time axis T and the periods P.
%
%   WAVESPECPLOT(T,X,P,W,R) optionally plots ABS(W).^R instead of W.
%  
%   WAVESPECPLOT(T,X,P,W,R,CI) makes a filled contour plot with contour
%   intervals CI.  If CI is not input, then the spectrum is plotted using 
%   PCOLOR, which is faster to render but slow to print. For making final 
%   figures, it is better to use the contouring option.
%
%   WAVESPECPLOT(T,X,P,W1,W2,...WN,...) makes an N+1 component plot, with
%   W1 in the second subplot, W2 in the third, etc.
%
%   After plotting, all subplots are packed together using PACKFIG.
%   H=WAVESPECPLOT(...) returns the handles to the subplots. 
%
%   [H,HL]=WAVESPECPLOT(...) also returns the handles to the lines plotted
%   in the first panel.
%   
%   If X is complex-valued, both the real and imaginary parts are plotted 
%   in the uppermost subplot.
%
%   Usage: h=wavespecplot(t,x,p,w);
%          h=wavespecplot(t,x,p,w,r);
%          [h,hl]=wavespecplot(t,x,p,w,r,ci);
% 
%   'wavespecplot --f' generates a sample figure.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1},'--f')
   type makefigs_wavespecplot;
   makefigs_wavespecplot
   return
end

t=varargin{1};
u=varargin{2};
p=varargin{3};
t=t(:);
p=p(:);

nspec=nargin-3;

n=1;
ci=[];

if length(varargin{end})==1
  if length(varargin{end-1})==1
     ci=varargin{end};
     n=varargin{end-1};
     nspec=nspec-2;
  else
     n=varargin{end};
     nspec=nspec-1;
  end
else
  if isvector(varargin{end})
     ci=varargin{end};
     n=varargin{end-1};
     nspec=nspec-2;
  end     
end

a=min(t);
b=max(t);

if isreal(u)
  c=maxmax(abs(u));
else
  c=max([maxmax(abs(real(u))) maxmax(abs(imag(u)))]);
end
c=c*1.1;

subplot(nspec+1,1,1)
if ~isreal(u)
  [hl1,hl2]=uvplot(t,u);
  hl=[hl1;hl2];
else
  hl=plot(t,u);
end
xlim([a,b])
ylim([-c c])
grid off
%hlines(0,'k:')

for i=1:nspec
  subplot(nspec+1,1,i+1)
  if ~isreal(varargin{3+i}) || n~=1
        varargin{3+i}=abs(varargin{3+i});
  end
  T=(varargin{3+i}').^n;
  if isempty(ci)
    pcolor(t,p,T),shading interp
  else 
    contourf(oprod(1+0*p,t),oprod(p,1+0*t),T,ci),hold on
    caxis([min(ci) max(ci)])
  end
  xlim([a,b])
  if p(1)<p(end)
      flipy
  end
  hold on
  set(gca,'tickdir','out')
  ylog
end
h=packfig(nspec+1,1,'rows');


