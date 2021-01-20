function[varargout]=gulf4plot(varargin)
%GULF4PLOT  A four-panel circulation plot for the Gulf of Mexico.
%
%   GULF4PLOT(LAT,LON,CV) makes a four-panel quiver plot in the Gulf of
%   Mexico with color shading being the speed of the mean flow.  CV is the 
%   complex-valued velocity U+iV of size LENGTH(LAT) x LENGTH(LON) x 4.
%
%   The plot has 2 rows and 2 columns, with a single colorbar centered 
%   beneath the four plots.  
%
%   GULF4PLOT(LAT,LON,CV,X) instead uses X, of the same size as CV, for the
%   color shading.
%
%   GULF4PLOT(...,CAX) uses CAX, a length 2 array, for the color axis.
%
%   GULF4PLOT(...,CAX,CLABEL,LABELS), where LABELS is a cell array of four
%   strings, prints LABELS{i} after the letter label of the ith subplot. 
%
%   LAT, LON, CV, and X can also be cell arrays of four elements.  This is
%   useful if the grid size is not the same in all four plots.
%
%   Usage: h=gulf4plot(lat,lon,cv);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details
 
% if strcmp(varargin{1}, '--t')
%     gulf4plot_test,return
% end

tweakmap_thin=['plot([-84.5+1i*21.5 -84.5+1i*22],''k--'',''linewidth'',0.5);hold on;'...
    'plot([-86.9+1i*21.5 -84.5+1i*21.5],''k--'',''linewidth'',0.5);' setstr(10) ...
    'jfig axis|[-98.35 -81.5 18 30.75] latratio|25 topoplot|gray '...
    'eval|topoplot([],[],-5:1:-1,''0.5E'') eval|topoplot([],[],-1/2,''0.5k'') eval|topoplot([],[],-0.005,''0.5E'') ' ...
    'ticksout portrait fontsize|[11 10 10 10] '];

cax=[];
clabel=[];
labels={'','','',''};
alpha=1/80;

if length(varargin{1})==1
    alpha=varargin{1};
    varargin=varargin(2:end);
end

lat=varargin{1};
lon=varargin{2};
cv=varargin{3};
x=varargin{4};
if length(varargin)>4
    cax=varargin{5};
end
if length(varargin)>5
    clabel=varargin{6};
end
if length(varargin)>6
    labels=varargin{7};
end

if isempty(x)
    if ~iscell(cv)
        x=abs(cv);
    else
        for i=1:length(cv)
            x{i}=abs(cv{i});
        end
    end
end

if ~iscell(x)
    xcell{1}=x(:,:,1);
    xcell{2}=x(:,:,2);
    xcell{3}=x(:,:,3);
    xcell{4}=x(:,:,4);
else
    xcell=x;
end

if ~iscell(cv)&&~isempty(cv)
    cvcell{1}=cv(:,:,1);
    cvcell{2}=cv(:,:,2);
    cvcell{3}=cv(:,:,3);
    cvcell{4}=cv(:,:,4);
else
    cvcell=cv;
end

if ~iscell(lat)
    latcell{1}=lat;
    latcell{2}=lat;
    latcell{3}=lat;
    latcell{4}=lat;
else
    latcell=lat;
end

if ~iscell(lon)
    loncell{1}=lon;
    loncell{2}=lon;
    loncell{3}=lon;
    loncell{4}=lon;
else
    loncell=lon;
end


figure
subplot(3,2,1)
jpcolor(loncell{1},latcell{1},xcell{1}),
hold on,eval(tweakmap_thin),
ar=get(gca,'dataaspectratio');
[long,latg]=meshgrid(loncell{1},latcell{1});
if ~isempty(cvcell)
    quiver([long(:);-90.2],[latg(:);20.5],alpha*[vcolon(real(cvcell{1}));30]*ar(1),alpha*[vcolon(imag(cvcell{1}));0],0,'k');hold on
end
text(-98,30.155,['(a) ' labels{1}])
text(-90.2,19.95,'30 cm/s')
%--------------------------------------------------------------------------
for i=2:4   
    subplot(3,2,i)
    jpcolor(loncell{i},latcell{i},xcell{i}),
    hold on,eval(tweakmap_thin),
    ar=get(gca,'dataaspectratio');
    [long,latg]=meshgrid(loncell{i},latcell{i});
    if ~isempty(cvcell)
        h=quiver(long,latg,alpha*real(cvcell{i})*ar(1),alpha*imag(cvcell{i}),0,'k');hold on
    end
    text(-98,30.155,['(' setstr(97+i-1) ') ' labels{i}])
end
%--------------------------------------------------------------------------
if ~isempty(cax)
    for i=1:4,subplot(3,2,i),caxis(cax),end
end
h1=packfig(3,2,'both');
delete(h1(5:6))
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
axes(h1(3)),set(gca,'xticklabelmode','auto')
axes(h1(4)),set(gca,'xticklabelmode','auto')

axes(h1(3))
hc=colorbar('South');
if ~isempty(clabel)
    hc.Label.String=clabel;
end
pos=hc.Position;
hc.Position=[0.32 0.34 0.4 pos(4)/2];
set(hc,'AxisLocation','out')

varargout{1}=h1;

 
function[]=gulf4plot_test
 
%reporttest('GULF4PLOT',aresame())
