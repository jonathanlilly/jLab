function[varargout]=ridgewalk(varargin)
% RIDGEWALK  Extract wavelet transform ridges, including bias estimates. 
%
%   [WR,IR,JR,FR]=RIDGEWALK(W,FS) where W is a wavelet transform matrix 
%   computed at frequecies FS, returns the wavelet ridges of transform W.
%
%   See WAVETRANS for details on the transform matrix W and frequencies FS.
%
%   RIDGEWALK returns the following quantities along ridges
%
%       WR     Wavelet transfrom value along the ridge
%       IR     Ridge indices into rows of W (time) 
%       JR     Ridge indices into columns of W (scale) 
%       FR     Instantaneous frequency along the ridge
%
%   All output variables are column vectors with the ridges appended one
%   atop the next, separated by a NaN.  Use COL2CELL(WR,IR,JR,FR) to 
%   convert these concatenated column vectors into cell arrays, or else
%   COL2MAT(WR,IR,JR,FR) to convert them into matrices.  
%
%   The minimum length of any ridge is two data points. 
%
%   The wavelet transform along the ridge, WR, estimates the analytic part 
%   of modulated oscillations present in original time series.  
%
%   RIDGEWALK(DT,...) uses a sample rate DT to compute the ridge frequency
%   FR.  The default value of DT is unity.  This does not affect the 
%   specification of FS, which is given in terms of a unit sample rate.
%   _______________________________________________________________________
%
%   Joint ridges
%
%   [W1R,W2R,...,WNR,IR,JR,FR]=RIDGEWALK(W1,W2,...,WN,FS) finds the joint
%   ridges of N transforms that all have the same size.  All output fields
%   remain column vectors.
%
%   In this case, there is only one set of ridges but N different transform
%   values. FR is then called the joint instantaneous frequency.
%
%   For details on joint ridges, see
%
%      Lilly and Olhede (2012), Analysis of Modulated Multivariate 
%           Oscillations. IEEE Trans. Sig. Proc., 60 (2), 600--612.
%   _______________________________________________________________________
%
%   Error estimate
%
%   [...,FR,ER]=RIDGEWALK(...,FS,P), where P=SQRT(BETA*GAMMA)
%   characterizes the generalized Morse wavelet used to form the wavelet 
%   transform, also returns an internal error estimate ER along the ridges.
%
%   This works for either univariate ridges or for the joint ridges.
%
%   ER measures the error with which the transforms estimate the analytic 
%   signals of modulated oscillations, arising from bias due to the 
%   modulation strength.  ER<<1 for signals that are accurately estimated. 
%
%   The expression for this quantity is given for a potentially 
%   multivariate signal by equation (62) of Lilly and Olhede (2012), and is 
%   based on a normalized version of the joint instantaneous curvature.   
%   _______________________________________________________________________  
%
%   Artifact removal
%   
%   RIDGEWALK has several features to minimize artifacts.
%
%   RIDGEWALK(...,FS,P,M) removes all ridges less than M*(2P/pi) periods in 
%   length. Since the number of periods in a generalized Morse wavelet is 
%   about 2P/pi, M gives the minimum number of wavelet lengths in a ridge.
%
%   To avoid spurious ridges due to the ridge analysis essentially seeing
%   the wavelet, one should definitely choose M>=1/2 and generally M>=1.
%   Experiments in noise show a big jump in ridge occurences below M=1/2. 
%
%   RIDGEWALK(...,FS,P,M,RHO) applies RIDGETRIM at level RHO, removing 
%   RHO*(P/pi) oscillations from the beginning and end of each ridge,
%   as these are generally contaminated by edge effects.  
%  
%   A choice of RHO=1 is recommended, or one wavelet half-width. 
%
%   RIDGETRIM is applied after the pruning set by M.  The shortest possible
%   ridge is then roughly (M-RHO)*(2P/pi).  Thus if the ridge trimming is 
%   applied, M will only have a net effect if it is greater than RHO.
%   _______________________________________________________________________  
%
%   Time-dependent frequency range
%
%   RIDGEWALK(...,[FMAX,FMIN]) specifies a maximum frequency and minumum
%   frequency FMAX and FMIN for the ridges.  These may be either scalars or
%   column vectors with the same number of rows as W.  Only ridge points 
%   between these two frequencies are used for the ridges.
%
%   FMAX and FMIN are both *radian* frequencies per unit time as specified
%   by DT, and thus have the same units as the ridge frequency FR. 
%   _______________________________________________________________________
%
%   Additional output arguments
%
%   [...,FR,ER,BR,CR]=RIDGEWALK(...,FS,P,...) additionally outputs the 
%   instantaneous bandwidth BR and curvature CR along ridges.  These may be
%   useful in error analysis, and are defined in equations (38) and (39) of 
%
%      Lilly and Olhede (2010), On the analytic wavelet transform. 
%           IEEE Trans. Info. Theory, 56 (8), 4135--4156.
%
%   For multivariate signals, BR and CR are *arrays* with one column per
%   input transform.  In this case, these quantities are defined by the
%   'deviation vectors' in equations (17) and (18) of Lilly and Olhede
%   (2012), divided by the squared norm of the wavelet transform. These 
%   reduce to the earlier definitions in the univariate case. 
%   _______________________________________________________________________  
%
%   See also WAVETRANS, RIDGEMAP, RIDGETRIM.
%
%   'ridgewalk --t' runs a test.
%   'ridgewalk --f' generates a sample figure.
%
%   Usage: [wr,ir,jr,fr,er]=ridgewalk(w,fs,P);
%          [wr,ir,jr,fr,er]=ridgewalk(w,fs,P,M);
%          [wr,ir,jr,wr,fr,er]=ridgewalk(dt,w,fs,P,M);
%          [wxr,wyr,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,P,M);
%          [wxr,wyr,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,P,M);
%          [wxr,wyr,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,P,M,rho);
%   _______________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2018 J.M. Lilly --- type 'help jlab_license' for details

%   Note that (62) of L&O (2012) reduces to (64) of L&O (2010)

%   For your information
%   _______________________________________________________________________
%
%   Interscale interpolation
%   
%   RIDGEWALK interpolates among discrete scale levels to yield more
%   accurate values for the ridge quantities WR and FR using a fast
%   quadratic interpolation.  
%   
%   See the low-level functions RIDGEINTERP and QUADINTERP for details.

%*************************************************************************   
%
%   Hidden options
%
%   Ridgewalk has several hidden options.  These are mostly use for 
%   testing purposes, so chances are you will not need to concern yourself
%   with them.  This documentation is provided for completeness. 
%
%   These options are the chaining parameter ALPHA, and the ridge type.
%   ___________________________________________________________________
%
%   ALPHA  -- Controls agressiveness of chaining across scales.
%
%   RIDGEWALK(...,'alpha',ALPHA) sets the chaining parameter ALPHA.
%
%   The chaining parameter ALPHA specifies the agressiveness with which
%   ridge points are chained across scales. The default value is one-half.
%
%   If desired, increase ALPHA to chain ridges more agressively across 
%   scales, or decrease ALPHA to supress chaining across scales.
%
%   ALPHA is defined as a normalized frequency difference
%
%         ALPHA  =  DOMEGA / OMEGA
%
%   where OMEGA is the transform frequency, and DOMEGA is the difference
%   between the frequency predicted for the next point based on the
%   transform at a "tail", and the actual frequency at prospective "heads".  
%
%   The chaining parameter is defined in such a way that it does not
%   need to be changed as time sampling or frequency sampling changes.
%   However, for strongly chirping signals or weakly chirping, noisy 
%   signals, better performance may perhaps be obtained by adjusting it.
%   ___________________________________________________________________
%
%   Ridge algorithms
%
%   RIDGEWALK(...,ALG)  determines the ridge algorithm.
%
%   Two different definitions may be used to locate the ridges.
%
%      ALG = 'phase'       Rate of transform change of phase definition
%      ALG = 'amplitude'   Maxima of transfom amplitude definition
% 
%   If ALG is not specified, 'amplitude' is used by default.
%
%   In practice, these usually do not differ much from one another.  An
%   examination of the difference between phase and amplitude ridges may be
%   found in 
%
%      Lilly and Olhede (2010).  On the analytic wavelet transform.
%
%   For reasons given therein, we prefer amplitude ridges.  This option
%   is therefore hidden, as the phase ridges are generally only used for 
%   testing purposes. Note our investigation disagreed with a comment in 
%   Mallat's book, which suggested that phase ridges should be superior.
%
%*************************************************************************


%   RIDGEWALK can be used to create a de-biased estimate of the signal,
%   following Lilly and Olhede (2010).  This estimate is given by
%
%       WR_DB= WR - (WR/2)*(BR^2+i*CR)*P^2
%
%   where WR, BR, and CR are defined above, where P is the wavelet 
%   time-frequency product, and i=SQRT(-1).  For the generalized Morse 
%   wavelets, P=SQRT(BETA*GAMMA).

%   This is now a hidden option.  No longer used by EDDYRIDGES
%   _______________________________________________________________________
%
%   Masked-out regions
%
%   RIDGEWALK permits the use to explicitly specify time-frequency regions 
%   which should be excluded from the ridge analyis.
%
%   RIDGEWALK(...,'mask',BOOL), where BOOL is a boolean array of the same
%   size as W, then those points for which BOOL is false will be excluded 
%   from the ridge analysis. In addition, ridges are not permitted to cross
%   such regions, to prevent spurious chaining between distant frequencies.
%
%   This functionality is useful if we have ancillary information, such as
%   a local signal-to-noise estimate, that can help determine a priori
%   which time-frequency points appear to be statistically significant. 


if strcmpi(varargin{1}, '--t')
    ridgewalk_test,return
elseif strcmpi(varargin{1}, '--f')
    type makefigs_ridgewalk
    makefigs_ridgewalk;
    return
end

%Default values of parameters
alg='amp';
rho=0;   %Trimming parameter
alpha=1/4;   %Ridge chaining parameter
P=[];
L=0;
chi=0;     %This is no longer changeable
fmax=[];
fmin=[];
mask=[];  

for i=1:max(3,length(varargin))
   if ischar(varargin{end})
        if strcmpi('amp',varargin{end}(1:3))||strcmpi('pha',varargin{end}(1:3))
            alg=varargin{end};
        end
        varargin=varargin(1:end-1);
    elseif ~ischar(varargin{end})&&ischar(varargin{end-1})
        if strcmpi('alp',varargin{end-1}(1:3))
            alpha=varargin{end};
        elseif strcmpi('mask',varargin{end-1}(1:3))
            mask=varargin{end};
        end
        varargin=varargin(1:end-2);
    end
end

if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else
    dt=1;
end

if size(varargin{end},2)==2
     fmax=varargin{end}(:,1).*dt;
     fmin=varargin{end}(:,2).*dt;
     varargin=varargin(1:end-1);
end

%Find all the W's
n=1;
while size(varargin{n+1})==size(varargin{1})
    n=n+1;
end
args=varargin(n+1:end);
varargin=varargin(1:n);

fs=args{1};
if length(args)==2
    P=args{2};
elseif length(args)==3
    P=args{2};
    L=args{3}*2*P./pi;
elseif length(args)==4
    P=args{2};
    L=args{3}*2*P./pi;
    rho=args{4};
end

[ir,jr,wr,fr,er,br,cr]=ridgewalk_one(varargin,P,L,alg,alpha,chi,dt,fmin,fmax,fs,mask);

%vsize(P,rho,fr,ir,jr,wr,er,br,cr)
if rho>0
     [fr,ir,jr,wr,er,br,cr]=ridgetrim(dt,P,rho,fr,ir,jr,wr,er,br,cr);
     disp(['RIDGEWALK trimming to ' int2str(length(ir)) ' ridge points.'])
end

%ir,jr,wr,fr,br,er
for n=1:length(varargin)
    if ~isempty(wr)
        varargout{n}=wr(:,n);
    else
        varargout{n}=[];
    end
end

varargout{n+1}=ir;
varargout{n+2}=jr;
varargout{n+3}=fr;
varargout{n+4}=er;
varargout{n+5}=br;
varargout{n+6}=cr;

disp('RIDGEWALK finished.')

function[ir,jr,wr,fr,er,br,cr]=ridgewalk_one(winput,P,L,alg,alpha,chi,dt,fmin,fmax,fs,mask)

[ir,jr,wr,fr,er,br,cr]=vempty;

if length(winput)>1
    w=zeros([size(winput{1},1) size(winput{1},2) length(winput)]);
    for i=1:length(winput)
        w(:,:,i,:)=winput{i};
    end
else
    w=winput{1};
end

if size(w,3)
    disp(['RIDGEWALK detecting a set of ' int2str(size(w,3)) ' transforms.'])
end

[bool,rq,wjoint,omjoint]=isridgepoint(w,fs,chi,alg,fmin,fmax,mask);

%figure,jpcolor(omjoint'),
%figure,jpcolor(diff(omjoint)')

disp('RIDGEWALK chaining ridges.')
[id,ir,jr,wr,fr]=ridgechains(fs,L,bool,wjoint,omjoint,alpha,mask,rq);
%figure,plot(diff(fr))

for i=1:nargout
    varargout{i}=[];
end

%/*************************************************************************
%Break ridges, interpolate, and compute bias parameters 
if ~isempty(id)
    [id,ir,jr,fr]=colbreaks(id,ir,jr,fr);

    %Account for sample rate
    fr=fr./dt;
    
    %Deviation vectors 
    w1=vdiff(w,1,'endpoint')./dt;
    w2=vdiff(w1,1,'nans')./dt;
    
    %Fix for the second derivative at the first and last point
    w2(1,:,:)=w2(2,:,:);
    w2(end,:,:)=w2(end-1,:,:);
    
    %Interpolate quantities within ridge for better estimation
    %Note fr has already been interpolated inside RIDGECHAINS
    [wr,w1r,w2r]=ridgeinterp(rq,ir,jr,w,w1,w2);
    
    %Multivariate formulation of the bias parameters
    frmat=vrep(fr,size(wr,2),2);
    wrmag=vrep(sqrt(sum(squared(wr),2)),size(wr,2),2);
    br=(w1r-sqrt(-1)*frmat.*wr)./wrmag;                     %See (17) of L012
    cr=(w2r-2*sqrt(-1)*frmat.*w1r-frmat.^2.*wr)./wrmag;     %See (18) of L012
    if ~isempty(P)
        er=frac(1,2)*squared(frac(P,fr)).*sqrt(sum(squared(cr),2));  %See (62) of LO12
    end
end
%\*************************************************************************


% The above reduces to this for a single time series
% [a,om,up,curv]=instmom(dt,w,'endpoint');
% [wr,fr,br,cr]=ridgeinterp(fs,rq,ir,jr,w,om,up,curv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[structout]=splitcells(x)
for i=1:length(x)
    for j=1:size(x{1},2)
       structout{j}{i,1}=x{i}(:,j);
    end
end

function[]=ridgewalk_test

load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
P=sqrt(2*4);
[wx,wy]=wavetrans(real(cx),imag(cx),{2,4,fs,'bandpass'},'mirror');
[wp,wn]=vectmult(tmat,wx,wy);

%Form ridges of component time series
[wr,ir,jr,fr,er]=ridgewalk(dt,wn,fs,P,1,'phase'); 
[wr2,ir2,jr2,fr2,er2]=ridgewalk(dt,wn,fs,P,1,'amplitude');

err=vsum(abs(wr-wr2).^2,1)./vsum(abs(wr).^2,1);
reporttest('RIDGEWALK phase and amplitude signal estimation error for NPG-06',err<1e-3)

[wpr,wnr,ir,jr,fr,er]=ridgewalk(dt,wp,wn,fs,P,3);   
[wxr,wyr,ir2,jr2,fr2,er2]=ridgewalk(dt,wx,wy,fs,P,3);   
[wpr2,wnr2]=vectmult(tmat,wxr,wyr);

err=vmean((abs([wpr wnr]-[wpr2 wnr2])./sqrt(abs([wpr wnr]).^2+abs([wpr2 wnr2])).^2),1);
reporttest('RIDGEWALK XY vs. PN invariance joint ridge for NPG-06',allall(err<1e-10))

%[wx,wy]=wavetrans(real(cx)+5*randn(size(cx)),imag(cx)+5*randn(size(cx)),{2,4,fs,'bandpass'},'mirror');
%[wxr,wyr,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,P,0); 
%[kappa,lambda,theta,phi]=ellparams(wxr,wyr);
%[x2,y2]=ellsig(kappa,lambda,theta,phi);
%aresame(wxr,x2,1e-12),aresame(wyr,y2,1e-12)
%These are identical, so I don't need to carry around the analytic parts
