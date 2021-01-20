function[ra,rb]=mconf(K,gamma,str)
%MCONF  Confidence intervals for the multitaper spectral estimate.
%
%   [RA,RB]=MCONF(K,GAMMA) returns the level GAMMA confidence interval 
%   for the direct multitaper power spectral estimate of an average over K
%   eigenspectra.  RA and RB lower and upper ratios to the true value.
%
%   More specifically, if S0 is the true value of the spectrum, and if S is
%   the spectral estimate, then RA and RB satisfy
%
%         Probability that RA < S/S0 < RB = GAMMA
%
%   with RA and RB defined to be symmetrically placed around 1. 
%
%   For example, GAMMA=0.95 computes the 95% confidence interval.
%
%   For the default settings of the multitaper spectrum implemented by
%   MSPEC, K is equal to 2P-1 where P is the time-bandwidth product.
%
%   The estimated confidence intervals can then be plotted as 
%
%         plot(f,S),hold on,plot(f,S*ra),plot(f,S*rb)
%
%   where S is the spectral estimated computed by MSPEC, while is the 
%   Fourier frequencies. 
%
%   MCONF relies upon the fact that S/S0 is approximately distributed
%   as 1/(2K) times a chi-squared distribution with 2K degrees of freedom.
%   _____________________________________________________________________
%  
%   Logarithmic confidence intervals  
%
%   When plotting the logarithm of the spectrum, one should use a different
%   set of confidence intervals.
%
%   MCONF(K,GAMMA,'log10') returns the confidence intervals for the base-10
%   logarithm of the multitaper spectral estimate.  In this case, RA and RB
%   are defined such that
% 
%         Probability that RA < LOG10(S)/LOG10(S0) < RB = GAMMA
%
%   and the confidence intervals can be plotted using either
%
%         plot(f,S),hold on,plot(f,10.^ra*S),plot(f,10.^rb*S)
%         set(gca,'yscale','log')
% 
%   or else 
%
%         plot(f,log10(S)),hold on,plot(f,ra+log10(S)),plot(f,rb+log10(S)).
%
%   MCONF(K,GAMMA,'natural') similarly returns the confidence intervals for
%   the natural logarithm of the multitaper spectral estimate.
%   _____________________________________________________________________
%
%   See also MSPEC, CHISQUARED. 
%
%   'mconf --t' runs a test.
%   'mconf --f' generates a sample figure.
%
%   Usage: [ra,rb]=mconf(K,gamma);
%          [ra,rb]=mconf(K,gamma,'log10');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details

   
%   For K greater than about 7, the standard and logarithmic confidence
%   intervals are quite similar when viewed on a logarithmic plot.  For 
%   smaller values of K

if strcmp(K, '--t')
    mconf_test,return
elseif strcmp(K, '--f')
    makefigs_mconf,return
end

if nargin==2
    str='linear';
end

% x=[0:0.0001:2]'*2*K;
% fx=chisquared(x,2*K);
% 
% [~,ii]=min(abs(x-2*K)); 
% 
% sumfx=cumsum(fx).*(x(2)-x(1));
% sumfx=abs(sumfx-sumfx(ii));

if strcmpi(str(1:3),'lin')
    dx=0.0001;
    %Compute pdf symmetrically about unity
    x1=[1-dx/2:-dx:-1]'*2*K;%from one to minus one
    x2=[1+dx/2:dx:3]'*2*K;  %from 1 to three 
    fx=chisquared(x1,2*K)*2*K+chisquared(x2,2*K)*2*K;
    sumfx=cumsum(fx).*dx;
    
    ii=find(sumfx>=gamma,1,'first');
    ra=x1(ii)/2/K;
    rb=x2(ii)/2/K;
elseif strcmpi(str(1:3),'nat')
    dx=0.0001;
    xo=log(2)+log(1/2/K)+psi(K); %see pav15-arxiv
    logx1=[xo-dx/2:-dx:xo-2]';
    logx2=[xo+dx/2:dx:xo+2]';
    x1=(exp(1).^logx1);
    x2=(exp(1).^logx2);
    fx=x1.*chisquared(2*K*x1,2*K)*2*K+x2.*chisquared(2*K*x2,2*K)*2*K;
    sumfx=cumsum(fx).*dx;
    ii=find(sumfx>=gamma,1,'first');
    ra=logx1(ii);%-log(2*K);
    rb=logx2(ii);%-log(2*K);
elseif strcmpi(str(1:3),'log')
    dx=0.0001;
    xo=log(2)+log(1/2/K)+psi(K);%see pav15-arxiv
    c=log(10);
    xo=xo./c;%change of base rule
    logx1=[xo-dx/2:-dx:xo-2]';
    logx2=[xo+dx/2:dx:xo+2]';
    x1=(10.^logx1);
    x2=(10.^logx2);
    fx=c*(x1.*chisquared(2*K*x1,2*K)*2*K+x2.*chisquared(2*K*x2,2*K)*2*K);

    sumfx=cumsum(fx).*dx;
    
    ii=find(sumfx>=gamma,1,'first');
    ra=logx1(ii);%-log(2*K);
    rb=logx2(ii);%-log(2*K);
end


function[]=mconf_test

rng(0);
P=4;
z=randn(2^20,1);
[psio,lambda]=sleptap(2^20,P);
[f,spp]=mspec(z,psio);   
[n1,x1]=hist(spp,1000);

K=2*P-1;

%fx=chisquared(x1*2*K,2*K)*2*K
%figure,plot(x1,n1./sum(n1)/(x1(2)-x1(1))),hold on, plot(x1,fx) 
%fx2=frac(K^K,factorial(K-1)).*x1.^(K-1).*exp(-K.*x1);plot(x1,fx2,'g.')

%test of confidence intervals for linear values of spectrum
gamma=length(find(spp>0.5&spp<1.5))./length(spp);  %This is a probability
[ra,rb]=mconf(K,gamma);
reporttest('MCONF linear with P=4',aresame(ra,0.5,1e-3)&&aresame(rb,1.5,1e-3))

%about the same gamma level, but for the natural log of the spectrum 
xo=log(2)+log(1/2/K)+psi(K);%psi here is the digamma function 
gamma=length(find(log(spp)>xo-.53&log(spp)<xo+.53))./length(spp);  %This is a probability
[ra,rb]=mconf(K,gamma,'natural');
reporttest('MCONF natural log with P=4',aresame(ra,xo-.53,1e-2)&&aresame(rb,xo+.53,1e-2))

%about the same gamma level, but for the natural log of the spectrum 
xo=(log(2)+log(1/2/K)+psi(K))/log(10);%psi here is the digamma function 
gamma=length(find(log10(spp)>xo-.23&log10(spp)<xo+.23))./length(spp);  %This is a probability
[ra,rb]=mconf(K,gamma,'log10');
reporttest('MCONF log10 with P=4',aresame(ra,xo-.23,1e-2)&&aresame(rb,xo+.23,1e-2))

function[]=mconf_figure



P=4;
psi=sleptap(length(cv),P);
[f,spp,snn,spn]=mspec(cv,conj(cv),psi); 
figure
h=plot(f,spp);linestyle -h h 2T
[ra,rb]=mconf(2*P-1,0.95);
hold on,plot(f,spp.*ra,'k'),plot(f,spp.*rb,'k')
[ra,rb]=mconf(2*P-1,0.95,'log10');
%hold on,plot(f,10.^(log10(spp)+ra),'r'),plot(f,10.^(log10(spp)+rb),'r')
hold on,plot(f,10.^(ra).*spp,'r'),plot(f,10.^(rb).*spp,'r')
ylog

figure
h=plot(f,log10(spp));linestyle -h h 2T
[ra,rb]=mconf(2*P-1,0.95,'log10');
hold on,plot(f,log10(spp)+ra,'r'),plot(f,log10(spp)+rb,'r')

% 
% 
% figure,plot(f,log10(spp)),hold on
% hold on,plot(f,log10(spp)+ra,'r'),plot(f,log10(spp)+rb,'r')
% 
% 
% 
% h=twospecplot(f,spp,snn);
% axes(h(1)),vlines(abs(corfreq(lat))), linestyle D b r k:
% ax=axis;axis([10^-2.95 ax(2) 10^-3 10^5]),xtick(10.^[-3 -2 -1 0])
% xlabel('Frequency (rad/hour)'),ylabel('Power Spectral Density')
% axes(h(2)),vlines(abs(corfreq(lat))), linestyle D b r k:
% ax=axis;axis([10^-2.95 ax(2) 10^-3 10^5]),xtick(10.^[-3 -2 -1 0])
% xlabel('Frequency (rad/hour)')



