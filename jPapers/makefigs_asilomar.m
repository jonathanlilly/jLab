function[varargout]=makefigs_asilomar(str)
%MAKEFIGS_ASILOMAR  Make figures for Lilly and Olhede (2009b).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M., and S. C. Olhede (2009). Wavelet ridge estimation of 
%      jointly modulated multivariate oscillations. IEEE 43rd Asilomar 
%      Conference on Signals, Systems, and Computers, 452--456. 
%
%   Usage: makefigs_asilomar
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

%/******************************************************************
load ebasnfloats
use ebasnfloats

len=cellength(lat);
index=find(len>200);
lato=30;
lono=-25;

vindex(id,num,lat,lon,p,t,index,1);


clear ga be fs p
for i=1:length(lat)
    if i==25||i==26
        %Special treatment for 36 (148000, now 25) and 37 (149000, now 26), 
        %       which are not getting ridges otherwise
        ga(i)=3;
        be(i)=8;
        fs{i}=morsespace(ga(i),be(i),{0.05,pi},2*pi/100,8);
    else
        ga(i)=3;
        be(i)=3;
        fs{i}=morsespace(ga(i),be(i),{0.05,2*pi/3},2*pi/100,8);
    end
    ppsi(i)=morseprops(ga(i),be(i));
end


%Compute wavelet transforms using generalized Morse wavelets
clear wx wy ir jr wxr wyr fxr fyr
for i=1:length(lat)
    cx{i}=fillbad(latlon2xy(lat{i},lon{i},lato,lono));   
    [wx{i},wy{i}]=wavetrans(real(cx{i}),imag(cx{i}),{ga(i),be(i),fs{i},'bandpass'},'mirror');
    dt=num{i}(2)-num{i}(1);
    [wxr{i},wyr{i},ir{i},jr{i},fr{i}]=ridgewalk(dt,wx{i},wy{i},fs{i},ppsi(i),3); 
end


clear xhatr yhatr xres lathat lonhat latres lonres kap lam the phi
clear fxhat fyhat cxhat cyhat fbar 

for i=1:length(lat)
    %Small overlap of two weak ridges for i=9 and i=23 
    [xhatr{i},yhatr{i},fbar{i}]=...
         ridgemap(length(cx{i}),wxr{i},wyr{i},fr{i},ir{i},'collapse');
    xhatmag=sqrt(squared(xhatr{i})+squared(yhatr{i}));
    [kap{i},lam{i},the{i},phi{i}]=ellparams(xhatr{i},yhatr{i});  
    xhat{i}=real(xhatr{i})+sqrt(-1)*real(yhatr{i});
    xres{i}=cx{i}-vswap(xhat{i},nan,0);
    [latres{i},lonres{i}]=xy2latlon(xres{i},lato,lono);
end

clear cv cvres cvhat
for i=1:length(lat)
    cv{i}=latlon2uv(num{i},lat{i},lon{i});
    cv{i}=fillbad(cv{i});
    cvres{i}=latlon2uv(num{i},latres{i},lonres{i});
    cvres{i}=fillbad(cvres{i});
    cvhat{i}=cv{i}-cvres{i};
end
jj=find(cellength(ir)>0);
%\******************************************************************

%/******************************************************************
figure
ax=[-30.5 -19.5 18 36.5];

subplot(131),h=cellplot(lon,lat);linestyle -h h D
hold on,h=cellplot(lon(jj),lat(jj));linestyle -h h k
h1=plot(lon{24},lat{24});linestyle -h h1 5C
h1=plot(lon{24},lat{24});linestyle -h h1 1k


ylabel('Latitude')
xlabel('West Longitude'),latratio(30),axis(ax)
title('Float Trajectories')

subplot(132)
xlabel('West Longitude'),title('Signal Ellipses')

for i=1:length(lat)
    if anyany(~isnan(kap{i}))
        dt=round(frac(2*pi,fbar{i}));
        clear index
        index(1)=dt(find(~isnan(dt),1,'first'));

        while index(end)<length(kap{i})-dt(find(~isnan(dt),1,'last'))
            index(end+1)=index(end)+dt(index(end));
        end
        index=nonnan(index);
        index=index(~isnan(kap{i}(index)));
        if ~isempty(index)
            ar=[frac(360,2*pi*radearth) frac(360,2*pi*radearth*cosd(30))];
            h=ellipseplot(2*kap{i},lam{i},the{i},lonres{i}+sqrt(-1)*latres{i},ar,'index',index);hold on
        end
    end
end
latratio(30),axis(ax)
linestyle k D

subplot(133),h=cellplot(lonres,latres);linestyle -h h D
hold on,h=cellplot(lonres(jj),latres(jj));linestyle -h h k
xlabel('West Longitude'),latratio(30),axis(ax),title('Residuals')
h1=plot(lonres{24},latres{24});linestyle -h h1 5C
h1=plot(lonres{24},latres{24});linestyle -h h1 1k

for i=1:3
    subplot(1,3,i)
    xtick([-30:2.5:-20])
    xtl=get(gca,'xticklabel');
    set(gca,'xticklabel',xtl(:,2:end))
    %fixlabels([-1,0])
    boxon
end
letterlabels(4)
packfig(1,3,'columns')


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 11 7])

if strcmpi(str,'print')
   print -depsc ebasn-trajectories.eps
end
%\******************************************************************


%/*************************************************    
figure
ci=(0:5:65);
ii=24;
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num{ii}-numo,cv{ii},2*pi./fs{ii},sqrt(abs(wx{ii}).^2+abs(wy{ii}).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Speed (cm/s)'),title('Bivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on
plot(num{ii}-numo,2*pi./fbar{ii},'w','linewidth',4)
plot(num{ii}-numo,2*pi./fbar{ii},'k','linewidth',2)

xlabel('Day of Year 1986'),ylabel('Period in Days')
set(gca,'ytick',2.^(2:.5:5.5))
set(gca,'yticklabel',[' 4';'  ';' 8';'  ';'16';'  ';'32';'  '])
inticks
text(-90,4.5,'(b)')


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 9 5])
if strcmpi(str,'print')
    print -deps ebasn-example.eps
end
%\*************************************************

%END of jlab_makefigs_asilomar
