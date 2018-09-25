function[varargout]=ridgetrim(varargin)
%RIDGETRIM  Trim edge effect regions from wavelet ridges.
%
%   Wavelet ridges, such as those computed by RIDGEWALK, are generally 
%   contaminated by edge effects at both ends.  In these locations, the
%   ridges give biased estimates of the signal properties. 
%   
%   FRO=RIDGETRIM(P,RHO,FR), where P=SQRT(BETA*GAMMA) and FR is a column
%   vector of frequencies from the ridge analysis, removes RHO*P/pi
%   periods from the beginning and end of each ridge.  Ridges containing
%   less than two points after this trimming are removed.
%
%   Note that P=SQRT(BETA*GAMMA) is characterizes the generalized Morse 
%   wavelet used in the wavelet transform.
%
%   P/pi is approximately the number of oscillations within the central
%   window of the wavelet, see eqn. (37) of Lilly and Olhede (2009).
%
%   The choice RHO=1 thus corresponds to trimming one wavelet half-width 
%   from end.  Numerical experiments show RHO=1 is a good value, but this
%   depends somewhat on the values of GAMMA and BETA chosen.
%
%   For joint ridges for multivariate signals FR should be a column vector 
%   of the joint instantaneous frequency.  This is formed from the array FR
%   output by RIDGEWALK by the weighted average FR=VMEAN(FR,1,ABS(WR).^2).
%
%   RIDGETRIM(DT,P,RHO,FR) specifies that FR has been computed with the
%   sample rate DT.
%
%   [FRO,X1O,X2O,...,XNO]=RIDGETRIM(P,RHO,FR,X1,X2,...,XN), where the XN
%   are any other ridge quantities, trims these in the same way as FR. The
%   XN are arrays having the same number of rows as FR.  
%
%   See also RIDGEWALK.
%   
%   'ridgetrim --f' generates a sample figure.
%
%   Usage: fr=ridgetrim(P,rho,fr);
%          [fr,wr,ir,jr]=ridgetrim(P,rho,fr,wr,ir,jr);
%          [fr,wr,ir,jr]=ridgetrim(dt,P,rho,fr,wr,ir,jr);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--f')
    type makefigs_ridgetrim
    makefigs_ridgetrim;
    return
elseif strcmp(varargin{1}, '--t')
   ridgetrim_test
   return
end
 
dt=1;
if length(varargin{3})==1
    dt=varargin{1};
    varargin=varargin(2:end);
end
P=varargin{1};
rho=varargin{2};
fr=varargin{3};
args=varargin(4:end);

[fr,args]=trimedges(dt,P,rho,fr,args);
varargout{1}=fr;
for i=1:length(args)
    varargout{i+1}=args{i};
end

function [fr,args]=trimedges(dt,P,rho,fr,args)

if ~isempty(fr)
    fr=col2cell(fr);
    
    for i=1:length(args)
        args{i}=col2cell(args{i});
    end
    
    for j=1:length(fr)        
        ar=cumsum(fr{j}./(2*pi).*dt,1);  %Ridge age in cycles
        br=ar(end)-ar;                   %Reverse ridge age

        a=find(ar>(rho*P/pi),1,'first');
        b=find(br>(rho*P/pi),1,'last');
        
        fr{j}=fr{j}(a:b);
        
        for i=1:length(args)
            args{i}{j}=args{i}{j}(a:b,:);
        end
    end
    
    index=find(cellength(fr)>1);
    fr=fr(index);
    fr=cell2col(fr);
    
    for i=1:length(args)
        args{i}=args{i}(index);
        args{i}=cell2col(args{i});
    end
    
    %fr=cell2col(cellprune(fr,'quiet')); 
    %for i=1:length(args)
    %    args{i}=cell2col(cellprune(args{i},'quiet'));
    %end
end

function []=ridgetrim_test
%With minimum length
dt=1;
rng(0);
N=10000;
z=randn(N,1)+1i*randn(N,1);
fs=2*pi./(logspace(log10(10),log10(100),50)');
beta=10;gamma=3;
[wx,wy]=wavetrans(real(z),imag(z),{gamma,beta,fs,'bandpass'},'mirror');
[wxr,wyr,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,sqrt(beta*gamma),2,1);
len=ridgelen(dt,fr);len=cellfirst(col2cell(len));
bool(1)=min(len)./(2*sqrt(beta*gamma)/pi)>1;
%min(len)./(2*sqrt(beta*gamma)/pi)

beta=3;gamma=3;
[wx,wy]=wavetrans(real(z),imag(z),{gamma,beta,fs,'bandpass'},'mirror');
[wxr,wyr,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,sqrt(beta*gamma),1.5,1);
len=ridgelen(dt,fr);len=cellfirst(col2cell(len));
bool(2)=min(len)./(2*sqrt(beta*gamma)/pi)>0.5;
%min(len)./(2*sqrt(beta*gamma)/pi)

reporttest('RIDGETRIM minimum ridge length is (M-RHO)*2*P/pi', bool)



