function[varargout]=makefigs_stokes(str)
%MAKEFIGS_STOKES  Make figures for Lilly et al. (2024).
%
%   This function makes all figures for the paper 
%
%   Lilly, Feske, Fox-Kemper, and Early (2024).  Integral theorems for the 
%       gradient of a vector field, with a fluid dynamical application.  
%       Proceedings of the Royal Society of London, Series A. 480 (2293): 
%       20230550, 1–30. https://doi.org/10.1098/rspa.2023.0550.
%
%   To run it, you will need the Crameri Perceptually Uniform Colormaps,
%   https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps.
%
%   Usage: makefigs_stokes
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2024 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

%This is where you want things to be printed to
dirname='/Users/lilly/Desktop/Dropbox/Projects/stokes/figures';

%This is when your NetCDF file is kept
ncdir='/Users/lilly/Desktop/Dropbox/NetCDF/';
filename='BetaEddyOne.nc';

set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultColorbarTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultTextInterpreter','latex')

%/************************************************************************
cmap=crameri('broc');
c1=cmap(1+40,:);

x=[-7:1:7]';y=x;
[xg,yg]=meshgrid(x,y);

%generate a bean curve
phi=pi*(0:1/256:1)';
z=(sin(phi).^1+cos(phi).^9).*exp(1i*phi);
z=8*(z-0.4*1i-.25)*exp(-1i*pi/8);

%exterior normal
dn=exp(-1i*pi/2)*vdiff(z,1,'periodic');

ii=1+(1:4:length(z)-1);

IJKL(:,:,1)=[1,0;0,1];
IJKL(:,:,2)=[0,-1;1,0];
IJKL(:,:,3)=[1,0;0,-1];
IJKL(:,:,4)=[0,1;1,0];

xvec=zeros(2,1,size(xg,1),size(xg,2));
xvec(1,1,:,:)=xg;
xvec(2,1,:,:)=yg;

u=zeros(2,1,size(xg,1),size(xg,2),4);
for j=1:4
    u(:,:,:,:,j)=pagemtimes(IJKL(:,:,j),xvec);
end

figure
for i=1:4
    for j=1:4
       subplot(4,4,j+(i-1)*4), axis equal,axis([-1 1 -1 1]*max(x)),
       %vlines(0,'k:'),hlines(0,'k:'), 
       noxlabels,noylabels,hold on,boxon,
       set(gca,'xtick','','ytick','')
       plot(z,'k','linewidth',2,'color',[1 1 1]*0.5)
       if i>0
           quiver(real(z(ii)),imag(z(ii)),real(dn(ii)),imag(dn(ii)),'color',[0.8500    0.3250    0.0980],'linewidth',0.5),
       end
       %set(gca,'xcolor',[1 1 1]*0.6,'ycolor',[1 1 1]*0.6)
       if i==1
            title(['$\mathbf{u}=\mathbf{' setstr(j+72) '}\mathbf{x}$     '],'interpreter','latex')
       end
       if j==1
           %if i==1
                %ylabel(['$\mathbf{u}$'],'interpreter','latex')%,'Rotation',0,'VerticalAlignment','middle')
           if i==2
                ylabel(['$\mathbf{' setstr(i+72) '}^T\mathbf{u}$       '],'interpreter','latex')%,'Rotation',0,'VerticalAlignment','middle')
           else
                ylabel(['$\mathbf{' setstr(i+72) '}\mathbf{u}$       '],'interpreter','latex')%,'Rotation',0,'VerticalAlignment','middle')
           end
       end
       if i==j
            set(gca,'linewidth',2,'xcolor','k','ycolor','k')
       end
    end
end

for i=1:4
    for j=1:4
       subplot(4,4,j+(i-1)*4),
       uij=pagemtimes(IJKL(:,:,i)',u(:,:,:,:,j));
       if i==1
           lc=[0 0 0];
       else
           lc=[1 1 1]*0.7;
       end
       quiver(x,y,squeeze(uij(1,1,:,:)),squeeze(uij(2,1,:,:)),'color',lc,'linewidth',0.5),
    end
end
packfig(4,4,'both')

set(gcf,'paperposition',[1 1 9 8.6])
fontsize 12 12 12 12
%str='print';h
if strcmp(str,'print')
   jprint(dirname,'stokes-illustration','epsc')
%   jprint(dirname,'stokes-illustration','-r400')
end
%\************************************************************************

%/*************************************************************************
%eddy plan view with triangle paths
lato=ncread([ncdir filename],'latitude');
x=ncread([ncdir filename],'x')/1000;
y=ncread([ncdir filename],'y')/1000;
u=permute(squeeze(ncread([ncdir filename],'u')),[2 1 3])*100;%cm/s
v=permute(squeeze(ncread([ncdir filename],'v')),[2 1 3])*100;%cm/s
cv=u+1i*v; clear u v
zeta=permute(squeeze(ncread([ncdir filename],'zeta')),[2 1 3]);
N=permute(squeeze(ncread([ncdir filename],'nu')),[2 1 3]);
S=permute(squeeze(ncread([ncdir filename],'sigma')),[2 1 3]);

fo=abs(corfreq(lato))/3600; %rad/s
xo=0.5*(y(144)+y(145));
yo=0.5*(y(80)+y(81));
edgecolor=0.6*[1 1 1];
edgewidth=1;

[z,lt,ut]=trianglepath(.1,0,50,8);
z=z*rot(pi/12)+15-1i*120;
for i=1:length(ut)
     ut{i}=ut{i}*rot(pi/12)+15-1i*120;
     lt{i}=lt{i}*rot(pi/12)+15-1i*120;
end
 
[zo,lto,uto]=trianglepath(.1,0,40,6);
zo=zo-110+1i*20;
for i=1:length(uto)
     uto{i}=uto{i}-110+1i*20;
     lto{i}=lto{i}-110+1i*20;
end
 
figure
subplot(3,1,1),contourf(x-xo,y-yo,zeta(:,:,300)./fo*100,200)
subplot(3,1,2),contourf(x-xo,y-yo,N(:,:,300)./fo*100,200)
subplot(3,1,3),contourf(x-xo,y-yo,S(:,:,300)./fo*100,200)
for i=1:3
    subplot(3,1,i)
    nocontours,hold on
    axis equal
    axis([-166 439 -185 120])
    caxis([-.099 .099]*100)
    plot(z,'w','linewidth',2)
    plot(zo,'w','linewidth',2)
    plot(z,'color',0.4*[1 1 1],'linewidth',1)
    plot(zo,'color',0.4*[1 1 1],'linewidth',1)
    %cmocean('diff')
    crameri('broc')
end

%--------------------------------------------------------------------------
%interpolate onto triangles and find moments

[zetac,sigmac,nuc]=vzeros(length(lt),2);
for i=1:length(lt)
    %upper triangle
    xc=real(ut{i});yc=imag(ut{i});
    zc=interp2(x-xo,y-yo,cv(:,:,300),xc,yc);
    [zetac(i,1),~,sigmac(i,1),nuc(i,1)]=curvemoments(xc,yc,zc);

    %lower triangle
    xc=real(lt{i});yc=imag(lt{i});
    zc=interp2(x-xo,y-yo,cv(:,:,300),xc,yc);
    [zetac(i,2),~,sigmac(i,2),nuc(i,2)]=curvemoments(xc,yc,zc);
end

zetac=zetac./fo*100;
nuc=nuc./fo*100;
sigmac=sigmac./fo*100;

%--------------------------------------------------------------------------
%plot average values in each triange 

subplot(3,1,1)
for i=1:length(lt)
    scatter(real(mean(ut{i})),imag(mean(ut{i})),100,zetac(i,1),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
    scatter(real(mean(lt{i})),imag(mean(lt{i})),100,zetac(i,2),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
end
subplot(3,1,2)
for i=1:length(lt)
    scatter(real(mean(ut{i})),imag(mean(ut{i})),100,nuc(i,1),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
    scatter(real(mean(lt{i})),imag(mean(lt{i})),100,nuc(i,2),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
end
subplot(3,1,3)
for i=1:length(lt)
    scatter(real(mean(ut{i})),imag(mean(ut{i})),100,sigmac(i,1),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
    scatter(real(mean(lt{i})),imag(mean(lt{i})),100,sigmac(i,2),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)  
end

%--------------------------------------------------------------------------
%interpolate onto triangles and find moments

[zetaco,sigmaco,nuco]=vzeros(length(lto),2);
for i=1:length(lto)
    %upper triangle
    xc=real(uto{i});yc=imag(uto{i});
    zc=interp2(x-xo,y-yo,cv(:,:,300),xc,yc);
    [zetaco(i,1),~,sigmaco(i,1),nuco(i,1)]=curvemoments(xc,yc,zc);

    %lower triangle
    xc=real(lto{i});yc=imag(lto{i});
    zc=interp2(x-xo,y-yo,cv(:,:,300),xc,yc);
    [zetaco(i,2),~,sigmaco(i,2),nuco(i,2)]=curvemoments(xc,yc,zc);
end

zetaco=zetaco./fo*100;
nuco=nuco./fo*100;
sigmaco=sigmaco./fo*100;

%--------------------------------------------------------------------------
%plot average values in each triange 

subplot(3,1,1)
for i=1:length(lto)
    scatter(real(mean(uto{i})),imag(mean(uto{i})),100,zetaco(i,1),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
    scatter(real(mean(lto{i})),imag(mean(lto{i})),100,zetaco(i,2),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
end
subplot(3,1,2)
for i=1:length(lto)
    scatter(real(mean(uto{i})),imag(mean(uto{i})),100,nuco(i,1),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
    scatter(real(mean(lto{i})),imag(mean(lto{i})),100,nuco(i,2),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
end
subplot(3,1,3)
for i=1:length(lto)
    scatter(real(mean(uto{i})),imag(mean(uto{i})),100,sigmaco(i,1),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
    scatter(real(mean(lto{i})),imag(mean(lto{i})),100,sigmaco(i,2),'filled',...
        'markeredgecolor',edgecolor,'linewidth',edgewidth)
end

%--------------------------------------------------------------------------
%colorbars and labels
for i=1:3
    subplot(3,1,i),
    xlabel('East-West Position (km)')    
    ylabel('North-South Position (km)')    
end

hax=packfig(3,1,'rows');
for i=1:3
    axes(hax(i))
    hc=colorbar('EastOutside');
    %pos=get(hc,'position');
    %set(hc,'position')
    text(-155,-170,['(' setstr(96+i) ')'])
    hc.Label.Interpreter = 'latex';
    hc.Ticks=[-0.75:.25:.75]*10;

    if i==1
        hc.Label.String = 'Vorticity $\zeta/f \times 100$';
    elseif i==2
        hc.Label.String = 'Normal Strain $\nu/f \times 100$';
    elseif i==3
        hc.Label.String = 'Shear Strain $\sigma/f \times 100$';
    end
end


set(gcf,'paperposition',[1 1 8 10])
fontsize 10 10 10 10 
%str='print';h
if strcmp(str,'print')
   jprint(dirname,'vortex-sampling','jpeg','-r400')
end
%\*************************************************************************


%/*************************************************************************
%alongtrack line plots

%--------------------------------------------------------------------------
%interpolate along tracks
cvc=interp2(x-xo,y-yo,cv(:,:,300),real(z),imag(z),'spline');
cvco=interp2(x-xo,y-yo,cv(:,:,300),real(zo),imag(zo),'spline');

Zc=interp2(x-xo,y-yo,zeta(:,:,300)./fo*100,real(z),imag(z));
Zco=interp2(x-xo,y-yo,zeta(:,:,300)./fo*100,real(zo),imag(zo));

Nc=interp2(x-xo,y-yo,N(:,:,300)./fo*100,real(z),imag(z));
Nco=interp2(x-xo,y-yo,N(:,:,300)./fo*100,real(zo),imag(zo));

Sc=interp2(x-xo,y-yo,S(:,:,300)./fo*100,real(z),imag(z));
Sco=interp2(x-xo,y-yo,S(:,:,300)./fo*100,real(zo),imag(zo));

%--------------------------------------------------------------------------
%find alongtrack shear component

dpar=vdiff(z,1);
upar=real(cvc.*conj(dpar))./abs(dpar);
uper=imag(cvc.*conj(dpar))./abs(dpar);

dparo=vdiff(zo,1);
uparo=real(cvco.*conj(dparo))./abs(dparo);
upero=imag(cvco.*conj(dparo))./abs(dparo);

%plot(uparo),hold on,plot(upero)

%find alongtrack shear 
%dudx=vdiff(real(cvco),1)./abs(vdiff(zo,1))/1000/100/fo*100;
dudx=vdiff(upar,1)./abs(vdiff(z,1))/1000/100/fo*100;
dvdx=vdiff(uper,1)./abs(vdiff(z,1))/1000/100/fo*100;
dudxo=vdiff(uparo,1)./abs(vdiff(zo,1))/1000/100/fo*100;
dvdxo=vdiff(upero,1)./abs(vdiff(zo,1))/1000/100/fo*100;

dudx(abs(dudx)>20)=nan;
dvdx(abs(dvdx)>20)=nan;
dudxo(abs(dudxo)>20)=nan;
dvdxo(abs(dvdxo)>20)=nan;

%--------------------------------------------------------------------------
%plotting

cmap=crameri('broc');
%c1=cmap(1+48,:);
%c2=cmap(end-48,:);
c1=cmap(1+40,:);
c2=cmap(end-40,:);


figure
clear hl
for i=1:3
    subplot(3,1,i)
    plot(real(zo),dudxo,'color',0.4*[1 1 1]),hold on,plot(real(zo),dvdxo,'--','color',0*[1 1 1])
    hl(3)=plot(real(z),dudx,'color',0.4*[1 1 1]);hl(4)=plot(real(z),dvdx,'--','color',0*[1 1 1]);
    if i==1
         axis([-166 439 -26 9.9])
         text(-155,-23,['(' setstr(96+i) ')'])
    else
         axis([-166 439 -13.5 13.5])
         text(-155,-11.5,['(' setstr(96+i) ')'])

    end       
    vlines(0,'D:'),hlines(0,'D:')
end

subplot(3,1,1),
h(1)=plot(real(zo),Zco,'color',c1,'linewidth',2);
h(2)=plot(real(z),Zc,'color',c2,'linewidth',2);
%h(2)=plot(real(z),Zc);
%linestyle(h,'2T- 2U-')
hl(1)=h(1);hl(2)=h(2);
subplot(3,1,2),
%h(1)=plot(real(zo),Nco);
%h(2)=plot(real(z),Nc);
h(1)=plot(real(zo),Nco,'color',c1,'linewidth',2);
h(2)=plot(real(z),Nc,'color',c2,'linewidth',2);
%linestyle(h,'2T- 2U-')
subplot(3,1,3)
%h(1)=plot(real(zo),Sco);
%h(2)=plot(real(z),Sc);
%linestyle(h,'2T- 2U-')
h(1)=plot(real(zo),Sco,'color',c1,'linewidth',2);
h(2)=plot(real(z),Sc,'color',c2,'linewidth',2);
%linestyle(h,'2T- 2U-')


subplot(3,1,1)
for i=1:length(lt)
    hl(5)=plot(mean(real(ut{i})),zetac(i,1),'wo','markerfacecolor','k','linewidth',1);hold on
end
for i=1:length(lto)
    plot(mean(real(uto{i})),zetaco(i,1),'wo','markerfacecolor','k','linewidth',1),hold on
end
subplot(3,1,2)
for i=1:length(lt)
    plot(mean(real(ut{i})),nuc(i,1),'wo','markerfacecolor','k','linewidth',1),hold on
end
for i=1:length(lto)
    plot(mean(real(uto{i})),nuco(i,1),'wo','markerfacecolor','k','linewidth',1),hold on
end
subplot(3,1,3)
for i=1:length(lt)
    plot(mean(real(ut{i})),sigmac(i,1),'wo','markerfacecolor','k','linewidth',1),hold on
end
for i=1:length(lto)
    plot(mean(real(uto{i})),sigmaco(i,1),'wo','markerfacecolor','k','linewidth',1),hold on
end

for i=1:3
    subplot(3,1,i),
    xlabel('East-West Position (km)')
    if i==1
        ylabel('Vorticity $\zeta/f \times 100$');
        hleg=legend(hl,'Ship track in eddy','Ship track in filament',...
            'Parallel shear','Perpendicular shear','Cell averages',...
            'location','best');
    elseif i==2
        ylabel('Normal Strain $\nu/f \times 100$');
    elseif i==3
        ylabel('Shear Strain $\sigma/f \times 100$');
    end
end
set(hleg,'Position',[0.65 0.71 0.1206 0.0740])

packfig(3,1,'rows')
set(gcf,'paperposition',[1 1 8 7])
fontsize 10 10 10 10 
if strcmp(str,'print')
   jprint(dirname,'vortex-profiles','jpeg','-r400')
end
%\*************************************************************************



