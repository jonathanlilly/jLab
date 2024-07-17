function[varargout]=makefigs_analytic(str)
%MAKEFIGS_ANALYTIC  Make figures for Lilly and Olhede (2010b).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M., and S. C. Olhede (2010). On the analytic wavelet 
%      transform. IEEE Transactions on Information Theory, 56 (8),
%      4135--4156.
%
%   Usage: makefigs_analytic
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2010--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

%/*************************************************
figure
ga=[3 3 3 3 nan 1/3.5 1 6 18 nan];be=[1.5 3 7 81 9./ga(5:end)];
[p,skew,kurt]=morseprops(ga,be);

t=1:5000;t=t-mean(t);
[psi,psif]=vzeros(5000,5,'nan');
for i=1:length(be)
    [psi(:,i),psif(:,i)]=morsewave(5000,1,ga(i),be(i),2*pi/40);
end

textstr{1}='(1.5,3)';
textstr{2}='(3,3)';
textstr{3}='(7,3)';
textstr{4}='(81,3)';
textstr{5}='';
textstr{6}='(31.5,0.29)';
textstr{7}='(9,1)';
textstr{8}='(1.5,6)';
textstr{9}='(0.5,18)';
textstr{10}='';

clear ha
for i=1:10
    ha(i)=subplot(2,5,i);
    if i~=5&&i~=10
        psi(:,i)=psi(:,i)./maxmax(abs(psi(:,i)));
        uvplot(t./p(i).*(2*pi/40),psi(:,i)),hold on,plot(t./p(i).*(2*pi/40),abs(psi(:,i))),linestyle k E-- 2k
        ylim([-1 1.05]),xlim([-4.5 4.5]),
        if i==3
            %text(2,1.5,'Examples of Generalized Morse Wavelets','fontsize',14)
            title('Examples of Generalized Morse Wavelets')
            %    pos=get(get(gca,'title'),'position');
            %    set(get(gca,'title'),'position',pos-[0 .01 0])
        end
        text(-4.2,-.9,textstr{i})
        fixlabels([0 -1]),hlines(0,'k:'),
        if i<4,hlines([-1 1.05],'k'),end
        if i==1,vlines(-4.5,'k'),end
        if i==3,vlines(4.5,'k'),end
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    elseif i==5
        f=(0:5000-1)./126;
        plot(f,2*psif(:,1:4)),ylim([0 2.02]),xlim([0 2.5]), linestyle E-- 2k k 2E--
        ytick(0:.5:2),fixlabels([-1 -1])
        boxoff,set(gca,'xtick',[],'ytick',[]),vlines(1,'k:'),
        %title('Frequency Domain')
        legend('a','b','c','d')
    elseif i==10
        plot(f,2*psif(:,6:10)),ylim([0 2.02]),xlim([0 2.5]), linestyle E-- 2k k 2E--
        ytick(0:.5:2),xlabel('Frequency / Peak Frequency'),fixlabels([-1 -1])
        boxoff,set(gca,'xtick',[],'ytick',[]),vlines(1,'k:')
        %title('Frequency Domain')
        legend('f','g','h','i')
    end
end

letterlabels(ha,1);
packfig(2,5)

orient landscape
fontsize jpofigure
set(gcf,'paperposition',[1 1 9 4.5])

if strcmpi(str,'print')
   %
   print -depsc analytic-morsies.eps
end
%\******************************


%/************************************************************
N=20;
ga1=[(1:1:21)];
be1=[(1:1:21)];
[ga,be]=meshgrid(ga1,be1);
[p,skew,kurt]=morseprops(ga,be);
dcell=morsederiv(N,ga,be);

x=zeros([size(ga) N]);
for n=1:length(dcell);
   if iseven(n)
       fact=n/2;
   else
       fact=(n-1)/2;
   end
   x(:,:,n)=dcell{n}./((p.^2).^fact)./factorial(n);
end

figure
for i=1:size(x,1)
    for j=1:size(x,2)
           if ga1(j)<6
                plot(2:size(x,3),squeeze(abs(x(i,j,2:end))),'ko','markersize',5,'markerfacecolor','k'),hold on,ylog
           else 
                plot(2:size(x,3),squeeze(abs(x(i,j,2:end))),'k.','color',[.6 .6 .6]),hold on,ylog
           end
           ylim([10^-6 5]),xlim([1.8 20.5]),xtick([2 3 4 5 10 15 20]),xlog
    end
end

for i=1:size(x,1)
    for j=1:size(x,2)  
           if ga1(j)==3
%                 h=plot(4:size(x,3),squeeze(abs(x(i,j,4:end)))); linestyle -h h 2w
                 h=plot(4:size(x,3),squeeze(abs(x(i,j,4:end)))); linestyle -h h F

           end
    end
end

hlines(1,'k:'),title('Morse Wavelet Decay with $\beta>1$, $1<\gamma< 6$')
xlabel('Derivative Number at Peak Frequency'),ylabel('Normalized Magnitude')

fontsize jpofigure
set(gcf,'paperposition',[1 1 3.5 3.5])

if strcmpi(str,'print')
   print -depsc analytic-decay.eps
end
%\********************************************************  



%/*************************************************************************
load ebasnfloats
use ebasnfloats
num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
vindex(num,lat,lon,1:549,1);

num=num-datenum(1986,1,1)+1;

%vindex(num,lat,lon,1:400,1);
cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));
cv=latlon2uv(num,lat,lon);

dt=(num(2)-num(1));
mlat=vmean(lat(:),1);

fmax=abs(corfreq(mlat))*frac(dt*24,2); %One cycle per 2 inertial periods = 1 cycle per 2.6 day
fmin=abs(corfreq(mlat))*frac(dt*24,40);%One cycle per 40 inertial periods = 1 cycle per 53 days
fs=morsespace(3,3,fmax,fmin,8);

ga=[3 3 3];
be=[1.5 3 7];
p=morseprops(ga, be);
c=squared(2./squared(p));

clear ir jr xr fr ar br cr xra fra bra cra wx
    
for i=1:3
    wx(:,:,i)=wavetrans(real(cx),{1,ga(i),be(i),fs,'bandpass'},'mirror');
    [xr{i},ir{i},jr{i},fr{i},~,br{i},cr{i}]=ridgewalk(dt,wx(:,:,i),fs,2*morseprops(ga(i),be(i)),2);  
    br{i}=real(br{i}./xr{i});
    cr{i}=cr{i}./xr{i};
    [xra(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
end

if 0
    x=[real(xra(:,1)) imag(xra(:,1))];
    for i=1:8
        x=[x;flipud(x)];
    end
    audiowrite('eddysound.wav',x,8192/2);
end
    
figure
M=4;N=3;
ci=(0:5:65)/2;
  
for i=1:size(xra,2)
    subplot(M,N,i)
    plot(num,vdiff([  real(cx-xra(:,i)) real(xra(:,i)) real(cx)],1))
    yoffset 25,
    linestyle  1.5k  G k
    axis([-95 num(end) -10 70 ])
    if i==1,ylabel('Velocity (cm/s)'),end
    title(['Ridge Analysis with $\gamma=$' int2str(ga(i)) ', $\beta=$' num2str(be(i))])
    
    subplot(M,N,i+N)
    contourf(num,2*pi./fs,abs(wx(:,:,i)'),ci),caxis([1 25]),colormap gray,flipmap,flipy,ylog,hold on,nocontours
    plot(num(ir{i}(~isnan(ir{i}))),2*pi./fs(jr{i}(~isnan(ir{i}))),'w','linewidth',4)
    plot(num(ir{i}(~isnan(ir{i}))),2*pi./fs(jr{i}(~isnan(ir{i}))),'k','linewidth',2)
    ytick([2 4 8 16 32])
    if i==1,ylabel('Period in Days'),end
    
    subplot(M,N,i+2*N)
    plot(num,fra(:,i),'k'),hold on,plot(num,bra(:,i),'k--')
    axis([-95  num(end) -.5 1.9]),ytick([0:.5:1.5])
    if i==1,ylabel('$\omega(t)$ $\&$ $\upsilon(t)$'),fixlabels([0 -1]),end
    hlines(0,'k:')
    
    subplot(M,N,i+3*N)
    uvplot(num,frac(cra(:,i),fra(:,i).^2)),hold on,ytick([-.2:.1:.2]),fixlabels([0 -1])
    plot(num,frac(bra(:,i),fra(:,i)).^2),%plot(num,[abs(frac(cra(:,i),fra(:,i).^2)) -abs(frac(cra(:,i),fra(:,i).^2))])
    linestyle k G 2k     
    axis([-95  num(end) -.225 .225]),hlines([-1 1]*c(i),'k:')
    if i==1,ylabel('$\rho_2(t)$ $\&$ $\rho_1^2(t)$'),end  
    
    xlabel('Day of Year 1986')
end

for i=1:M*N
    h(i)=subplot(M,N,i);%text(,['(' setstr(i+96) ')' ])
    axes(h(i))
    xtick([-100:100:500])
    set(gca,'xticklabelmode','auto')
end
letterlabels(h,1);

packfig(4,3)
set(gcf,'paperposition',[1 1 9 4.5])
fontsize jpofigure
if strcmpi(str,'print')
   print -depsc analytic-transforms.eps
end

% %Reconstructing mean and median estimates
% for i=1:3
%      wx(:,:,i)=wavetrans(real(xra(:,i)),{1,ga(i),be(i),fs,'bandpass'},'mirror');
%     [ir{i},jr{i},xr{i},fr{i},br{i},cr{i}]=ridgewalk(dt,wx(:,:,i),fs,2*morseprops(ga(i),be(i)));  
%     [xra2(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
% end
% 
% for i=1:3
%     vmean(squared(xra2(:,i)-xra(:,i))./squared(xra(:,i)),1)
%     vmedian(squared(xra2(:,i)-xra(:,i))./squared(xra(:,i)),1)
% end

%END of jlab_makefigs_analytic