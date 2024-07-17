function[varargout]=isomax(varargin)
% ISOMAX  Returns those transform maxima that are isolated from others.
%  
%   This function is part of 'element analysis' described in Lilly (2017), 
%   "Element analysis: a wavelet-based method for analyzing time-localized
%   events in noisy time series", available at www.jmlilly.net.
%  
%   BOOL=ISOMAX(SIZ,INDEX,WW,FF,GAMMA,BETA,MU,LAMBDA) returns an boolean
%   array is true for transform maxima that are isolated from other maxima.
% 
%   SIZ is the original size of the wavelet transform array, while the 
%   index into maxima locations INDEX, transform values WW, and transform
%   frequencies FF, are all as computed by TRANSMAX. 
%
%   ISOMAX returns a boolean array BOOL that is true for those maxima that 
%   are isolated at level LAMBDA.  The maxima are interpreted as being from  
%   the transform *of* a (MU,GAMMA) wavelet *with* a (BETA,GAMMA) wavelet.
%
%   Specifically, LAMBDA defines a contour around a transform maxima at 
%   which the transform modulus is expected to have decayed to a fraction
%   LAMBDA of its peak value.  A maxima is 'isolated' if it is larger in 
%   magnitude that all other points encompassed that region of influence.
%
%   Note that if the third or fourth element of SIZ is greater than one, 
%   this indicates a set of original time series organized as a 2D or 3D 
%   array respectively, according to the WAVETRANS convention.  In this 
%   case the isolation criterion is applied only to events from the same 
%   time series; the extra dimensions are simply looped over.  
%
%   [BOOL,Z]=ISOMAX(...) also returns a matrix Z with LENGTH(INDEX) columns
%   containing the regions of influence T+1i*F surrounding all of the input
%   maxima.  One may then use 'plot(z)' to plot these regions.   
%
%   ISOMAX(...,'parallel') optionally uses a PARFOR loop when computing the
%   regions of influence, which can be useful for large datasets.
%
%   See also MAXPROPS, TRANSMAX, TRANSMAXDIST, MAX2EDDY.
%
%   Usage: bool=isomax(siz,index,ww,ff,gamma,beta,mu,lambda);
%          [bool,z]=isomax(siz,index,ww,ff,gamma,beta,mu,lambda);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017 J.M. Lilly --- type 'help jlab_license' for details        


if strcmpi(varargin{1},'--f')
%   This doesn't work right now.  
%   'isomax --f' generates a sample figure.
%    isomax_fig;
    return
end

parstr='serial';
str='band';
for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
            parstr=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'ban')||strcmpi(varargin{end}(1:3),'ene')
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end
        
siz=varargin{1};
index=varargin{2};
ww=varargin{3};
ff=varargin{4};
ga=varargin{5};
be=varargin{6};
mu=varargin{7};
lambda=varargin{8};

[ii,jj,nn,kk]=ind2sub(siz,index);
   
bool=false(size(ii));

if nargout == 2
    z=zeros(1000,length(ii));
end

if length(siz)==2
    siz(3)=1;
    siz(4)=1;
elseif length(siz)==3
    siz(4)=1;
end

[tz,fz]=morseregion(lambda,ga,be,mu,1,str);
%figure,plot(tz,fz)

for n=1:siz(3)
    %n
    booln=(nn==n);
    for k=1:siz(4)
        index=find(booln&(kk==k));
        %length(find(index))
        if ~isempty(index)
            if nargout ==1
                bool(index)=isomax1(ii(index),ww(index),ff(index),ga,be,mu,lambda,tz,fz,parstr,str);
            elseif nargout ==2
                [bool(index),z(:,index)]=isomax1(ii(index),ww(index),ff(index),ga,be,mu,lambda,tz,fz,parstr,str);
            end
        end
    end
end


varargout{1}=bool;
if nargout==2
    varargout{2}=z;
end

function[bool,z]=isomax1(ii,ww,ff,ga,be,mu,lambda,tz,fz,parstr,str)

originalposition=[1:length(ww)]';
%Sort by magnitude
if ~isempty(ww)
    [temp,sorter]=sort(-abs(ww));
    %sorter
    vindex(ww,ii,ff,originalposition,sorter,1);
end

bool=[];
[c,rho,frho]=maxprops(1,1,ga,be,mu,str);
fact=1./frho;
%fact=frac(morsefreq(ga,be),morsefreq(ga,mu)).*frac(mu+1,be).^(1./ga);
%aresame(fact,1./frho,1e-10)

if isempty(ga)||isempty(ii)
    z=[];
else
    z=vrep(tz+1i*fz,length(ff),2);
    fmat=vrep(ff(:)',length(tz),1);
    z=real(z)./fmat.*fact+1i*imag(z).*fmat./fact;
%     if strcmpi(parstr(1:3),'par')
%         parfor i=1:length(ff)
%             z(:,i)=tz./ff(i).*fact+sqrt(-1)*fz.*ff(i)./fact;
%         end
%     else
%         for i=1:length(ff)
%             z(:,i)=tz./ff(i).*fact+sqrt(-1)*fz.*ff(i)./fact;
%         end
%     end
%     
%     
    %z2=morseregion(fs(jj),lambda,ga,be,mu);
    %aresame(z,z2,1e-6)
    %plot(z,'k'),hold on,plot(z2,'r--')
    %z=real(z)+sqrt(-1)*imag(z); % figure,plot(z)
    %tmat=vrep(t(ii)',size(z,1),1);
    tmat=vrep(ii',size(z,1),1);
    %vsize(tmat,t,ii,z)
    z=z+tmat;
    bool=true(size(ii));

    if strcmpi(parstr(1:3),'par')
        parfor i=2:length(ww)
            ti=ii([1:i-1]);
            fi=ff([1:i-1]);
            bool(i)=~anyany(inpolygon(ti,fi,real(z(:,i)),imag(z(:,i))));
        end
    else
        for i=2:length(ww)
            ti=ii([1:i-1]);
            fi=ff([1:i-1]);
            bool(i)=~anyany(inpolygon(ti,fi,real(z(:,i)),imag(z(:,i))));
        end
    end

    bool(originalposition)=bool;
    z(:,originalposition)=z;
    %vindex(ii,jj,ff,ww,bool,1);
    %vindex(z,bool,2);
    %if ~isempty(ii)
    %    [ii,sorter]=sort(ii);
    %    vindex(jj,ff,ww,sorter,1);
    %    vindex(z,sorter,2);
    %end
end

%if isempty(ii)
%    disp('ISOMAX found no isolated maxima with the specified properties.')
%end

function[]=isomax_fig





z1=morseregion(10,.5,ga,be,mu);
z2=morseregion(70,.5,ga,be,mu);

    
%L=frho/2/pi;

ga=2;be=3;mu=0;frho=2*pi/100;
rho=morsefreq(ga,mu)./frho;
psi=morsewave(1000,ga,mu,frho)*rho;
psi=psi./maxmax(abs(psi));

fmax=frac(morsefreq(ga,be),rho*frac(be,mu+1).^(1./ga));

%maxmax(abs(psi))
%morseafun(ga,mu)./(2*pi*ga).*gamma((mu+1)./ga)./rho

%morseafun(ga,mu)./(2*pi*ga).*gamma((mu+1)./ga)./ maxmax(abs(psi))

num=[1:length(psi)]';
fs=morsespace(ga,be,2*pi/30,2*pi/1000,8);
s=morsefreq(ga,be)./fs';
%w=wavetrans(psi,{1,ga,be,fs,'bandpass'},'mirror')*sqrt(2);
w=wavetrans(real(psi),{ga,be,fs});
figure,
jpcolor(num,fs,abs(w)');shading interp,hold on,ylog,
[ii,jj,ww,ff]=transmax(fs,w);
plot(ii,ff,'m.')
[ii,jj,ww,ff,ahat,fhat,z]=isomax(ii,jj,ww,ff,ga,be,mu,0.1);
plot(ii,ff,'*'),plot(z)


figure,
lambda=0.9;%lambda=1/2;lambda=0.1;
[ii,jj,ww,ff]=transmax(fs,w);
[ii,jj,ww,ff,ahat,fhat,z]=isomax(ii,jj,ww,ff,ga,be,mu,lambda);
contour(num,fs,abs(w)'./maxmax(abs(w)),[lambda lambda],'k');hold on,plot(z)



[ix,jx,wx,ahat,fhat,z]=isomax(w,0,num,fs,1/2,ga,be,mu);colormap lansey
%plot(num(ix),fs(jx),'m*')


phio=frac(morseafun(ga,be).*morseafun(ga,mu),2*pi*ga).*...
    gamma(frac(be+mu+1,ga)).*...
    frac((s.^be).*(rho.^(mu+1)),(s.^ga+rho.^ga).^((be+mu+1)./ga));

phio=frac(morseafun(ga,be).*morseafun(ga,mu),2*pi*ga).*...
    gamma(frac(be+mu+1,ga)).*...
    frac((s.^be).*(rho.^(mu+1)),(s.^ga+rho.^ga).^((be+mu+1)./ga));


phio=morseafun(ga,be).*frac(gamma(frac(be+mu+1,ga)),gamma(frac(mu+1,ga))).*...
    frac((s.^be).*(rho.^(mu+1)),(s.^ga+rho.^ga).^((be+mu+1)./ga));





%fact=frac(morseafun(ga,be).*morseafun(ga,mu),morseafun(ga,be+mu)).*
%frac((s.^be).*(rho.^mu),(s.^ga+rho.^mu).^((be+mu+1)./ga));
%psio=morsewave(1000,be
%fact=vrep(fact,

%w,chi,t,fs,lambda,gamma,beta,mu


load bravo94
use bravo94.rcm
numo=datenum(1994,1,1)-1;

ga=2;be=1;mu=0;
fs=morsespace(ga,be,2*pi/30,2*pi/1000,8);

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cv(:,2)),imag(cv(:,2)),{ga,be,fs});
w=sqrt(abs(wx).^2+abs(wy).^2);
figure,
jpcolor(num-numo,24*fs,w');shading interp,hold on,ylog,colormap lansey
[ii,jj,ww,ff]=transmax(fs,w);
plot(num(ii)-numo,ff*24,'m.')
[ii,jj,ww,ff,ahat,fhat,z]=isomax(ii,jj,ww,ff,ga,be,mu,0.1);
plot(num(ii)-numo,ff*24,'k*')



figure,
jpcolor(1:size(w,1),24*fs,w');shading interp,hold on,ylog,colormap lansey
[ii,jj,ww,ff]=transmax(fs,w);
plot(ii,ff*24,'m.')
[ii,jj,ww,ff,ahat,fhat,z]=isomax(ii,jj,ww,ff,ga,be,mu,1/2);
plot(ii,ff*24,'k*')
plot(real(z)+1i*imag(z)*24)

figure
jpcolor(1:size(w,1),fs,w');shading interp,hold on,ylog,colormap lansey
[ii,jj,ww,ff]=transmax(fs,w);
plot(ii,ff,'m.')
[ii,jj,ww,ff,ahat,fhat,z]=isomax(ii,jj,ww,ff,ga,be,mu,1/2);
plot(ii,ff,'k*')
plot(z)
 
 

plot(z,'w')
[ix,jx,wx,z]=isomax(w,10,num-numo,fs*24,0.8,ga,be,0);
plot(z,'m')



[ix,jx,wx]=isomax(w,10,fs,ga,be,2);




[ix,jx,wx]=isomax(w,10,fs,ga,be,2,'plot',num-numo,fscale,'c');caxis([1 5])


%set(gca,'dataaspectratio',[200 1 1])
plot(num(ix)-numo,fscale*fs(jx),'m*')
xlabel('Day of Year 1994');
ylabel('Frequency (cycles per day)')
title('Isolated Maxima at the Bravo mooring')


psi=morsewave(10000,ga,be,2*pi/1000);
num=[1:length(x)]';num=num-mean(num);
w=wavetrans(x,{1,ga,be,fs,'bandpass'},'mirror');
figure,
jpcolor(num,fscale.*fs,sqrt(abs(w))');shading interp,hold on,ylog,flipmap
[ix,jx,wx]=isomax(w,10,fs,ga,be,2,'plot',num,fscale,'c');colormap lansey
plot(num(ix),fscale*fs(jx),'m*')

figure,
x=morsewave(10000,ga,be,2*pi/100);
fs=morsespace(ga,be,2*pi/10,2*pi/5000,8);
num=[1:length(x)]';num=num-mean(num);
w=wavetrans(x,{1,ga,be,fs,'energy'},'mirror');
%contourf(num,fs,log(abs(w))',100);nocontours,hold on,ylog,
contourf(num,fs,(abs(w))',100);nocontours,hold on,ylog,
%[ix,jx,wx]=isomax(w,1e-6,fs,ga,be,2,'plot',num,1,'c');colormap lansey
%plot(num(ix),fs(jx),'m*')

%figure,
%jpcolor(sqrt(abs(w))');shading interp,hold on,ylog,flipmap
%plot(ix,jx,'m*')

figure,jpcolor(num,fs,(abs(w')));,hold on,ylog,
[t,f]=morseregion(0.6,ga,be,be,2*pi/100);
plot(t,f)
[t,f]=morseregion(1,ga,be,2*pi/100);

figure,jpcolor(num,fs,abs(w)'),shading interp,hold on
contour(num,fs,(abs(w))',20,'k');
figure,jpcolor(num,flipud(fs),flipud(abs(w)')),shading interp,hold on
contour(num,flipud(fs),flipud(abs(w)'),20,'k');

figure,jpcolor(abs(w)'),shading interp,hold on
contour((abs(w))',20);

