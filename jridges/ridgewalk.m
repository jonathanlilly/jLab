function[varargout]=ridgewalk(varargin)
% RIDGEWALK  Extract wavelet transform ridges, including bias estimates. 
%
%   [IR,JR,WR,FR]=RIDGEWALK(W,FS) where W is a wavelet transform matrix
%   at frequecies FS, such as that returned by MORSEWAVE, returns the 
%   wavelet ridges of transform W.
%
%   The columns of W correspond to different frequencies, specified by the
%   frequency array FS, at which the wavelet transform was performed.  Note
%   that FS assumes a unit sample rate.  
%
%   The frequencies FS are expected to be ordered from highest to lowest.  
%
%   RIDGEWALK returns the following quantities along ridges
%
%       IR     Ridge indices into rows of W (time) 
%       JR     Ridge indices into columns of W (scale) 
%       WR     Wavelet transfrom value along the ridge
%       FR     Transform frequency values in radian frequency
%
%   All output variables are cell arrays with one ridge per cell.  To plot
%   the ridges, use CELLPLOT(IR,JR).  CELL2COL(IR,JR,WR,FR) turns these
%   cell arrays into concatenated column vectors. 
%
%   RIDGEWALK(DT,W,FS,...) uses a sample rate DT to compute the ridge
%   frequency FR.  The default value of DT is unity.  This does not affect
%   the specification of FS, which is given in terms of a unit sample rate.
%   _______________________________________________________________________
%
%   Masked-out regions
%
%   RIDGEWALK permits the use to explicitly specify time-frequency regions 
%   which should be excluded from the ridge analyis.
%
%   RIDGEWALK(...,W,FS,BOOL), where BOOL is a boolean array of the same
%   size as W, then those points for which BOOL is false will be excluded 
%   from the ridge analysis. In addition, ridges are not permitted to cross
%   such regions, to prevent spurious chaining between distant frequencies.
%
%   This functionality is useful if we have ancillary information, such as
%   a local signal-to-noise estimate, that can help determine a priori
%   which time-frequency points appear to be statistically significant. 
%   _______________________________________________________________________
%
%   Additional options
%
%   RIDGEWALK has other contingencies for rejecting spurious ridge points.
%   These tend to occur on the flanks of interesting signals, and 
%   reflect the wavelet structure rather than the signal structure.
%
%   RIDGEWALK(...,{L,CHI}) specifies options for the ridge computation.
%
%        L  -- Removes all ridges of less than L periods in length
%      CHI  -- Removes all small amplitude ridge points having |W|<CHI
%
%   In general, L should be set in proportion to the number of oscillations
%   contained in the wavelet.  A recommended setting is L > 2*P/pi, with P 
%   described shortly. This criterion means that the ridge must be longer
%   than the number of oscillations in the central envelope of the wavelet.
%
%   Here P is a quantity that characterizes the number of oscillations in a
%   wavelet.  For the generalized Morse wavelets calculated by MORSEWAVE, 
%   P is given by P=SQRT(BETA*GAMMA), see Lilly and Olhede (2009).
%
%   The options cell may also include some additional parameters for hidden
%   options that are used during testing. These should not be required by
%   most users, but are documented in the function body for completeness.
%   _______________________________________________________________________  
%
%   Time-dependent frequency range
%
%   The ridge curves may be limited to a time-varying frequency range.
%
%   RIDGEWALK(DT,...,{FMAX,FMIN,L,CHI}) additionally specifies a maximum 
%   frequency FMAX and minumum frequency FMIN for the ridges.  Only ridge
%   points between these two frequencies are returned.
%
%   FMAX and FMIN are both *radian* frequencies per unit time as specified
%   by DT. DT is optional and its default value is unity.  Thus FMAX and 
%   FMIN are directly comparable to the ridge frequency FR. 
%   
%   Both FMAX and FMIN are either scalars, or arrays the same length as W. 
%   _______________________________________________________________________
%
%   Output of bias parameters
%
%   [IR,JR,WR,FR,BR,CR]=RIDGEWALK(...) optionally outputs two additional
%   quantities along the ridges.
%
%       BR     Instantaneous bandwidth  
%       CR     Instantaneous curvature  
%
%   When these 'bias parameters' BR and CR are small compared with the
%   frequency, i.e. BR/FR << 1 and CR/(FR^2) << 1, then the signal is 
%   accurately estimated, as discussed in 
%
%      Lilly and Olhede (2010), On the analytic wavelet transform. 
%           IEEE Trans. Info. Theory, 56 (8), 4135--4156.
%
%   For more details, see INSTMOM.
%   _______________________________________________________________________
%
%   Joint ridges
%
%   [IR,JR,W1R,W2R,...,WNR,F1R,F2R,...FNR]=RIDGEWALK(W1,W2,...,WN,FS) 
%   finds the joint ridges of N transforms that all have the same size.  
%
%   In this case, there is only one set of ridges but N different values.
%   All of the output fields are again cell arrays of the same size.
%
%   [IR,JR,WR,FR]=RIDGEWALK(W1,W2,...,WN,FS,'mat') specifies an alternate,
%   'matrix' format for joint ridges.  Now WR and FR are cell arrays, each 
%   element of which is a matrix with N columns.  The default 'column' 
%   format splits these N columns into N different output arguments.
%
%   For details on joint ridges, see
%
%      Lilly and Olhede (2012), Analysis of Modulated Multivariate 
%           Oscillations. IEEE Trans. Sig. Proc., 60 (2), 600--612.
%   _______________________________________________________________________
%   
%   Bias parameters for joint ridges
%
%   [IR,JR,WR,FR,BR,CR]=RIDGEWALK(W1,W2,...,WN,FS,'mat') for the case of 
%   joint ridges similarly outputs the two bias parameters along the 
%   ridges, here organized in matrix format.
%
%   BR and CR are normalized versions of (17) and (18) of Lilly and Olhede 
%   (2012). Both BR and CR have been normalized by dividing them by the 
%   modulus of the estimated analytic signal, SQRT(SUM(ABS(WR).^2,2)).
%
%   [...,B1R,B2R,...,BRN,C1R,C2N,...,CNR]=RIDGEWALK(W1,W2,...,WN,FS), the
%   default, uses column format instead.  In this case there will be 2+4*N
%   total output arguments: IR, JR, and N columns of WR, FR, BR, and CR.     
%
%   The bias associated with the estimated signal is small compared to the 
%   magnitude of the signal when the modulus of XCR is small.
%
%   For more details on the bias parameters for multivariate signals see 
%   Lilly and Olhede (2012).
%   _______________________________________________________________________
%
%   Interscale interpolation
%   
%   RIDGEWALK interpolates among discrete scale levels to yield more
%   accurate values for the ridge quantities WR and FR using a fast
%   quadratic interpolation.  
%   
%   See the low-level functions RIDGEINTERP and QUADINTERP for details.
%   _______________________________________________________________________
%
%   See also WAVETRANS, RIDGEMAP.
%
%   'ridgewalk --t' runs a test.
%   'ridgewalk --f' generates a sample figure.
%
%   Usage: [ir,jr,wr,fr]=ridgewalk(w,fs);
%          [ir,jr,wr,fr]=ridgewalk(w,fs,{L,CHI});
%          [ir,jr,wr,fr]=ridgewalk(dt,w,fs,{L,CHI});
%          [ir,jr,wxr,wyr,fxr,fyr]=ridgewalk(dt,wx,wy,fs,{L,CHI});
%          [ir,jr,wr,fr,br,cr]=ridgewalk(dt,wx,wy,fs,{L,CHI},'mat');
%          [ir,jr,wr,fr]=ridgewalk(dt,w,fs,bool,{L,CHI});
%   _______________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details

%
%   Possibly say more about bias last
%         (1/2) P^2 SQRT(SUM(ABS(X2R).^2,2))./ 

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
%   RIDGEWALK(...,{N,CHI,ALPHA}) sets the chaining parameter ALPHA.
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
%   RIDGEWALK(...,{...,ALG}) where ALG is a string, and is the last
%   quanity in the options cell, determines the ridge algorithm.
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


if strcmpi(varargin{1}, '--t')
    ridgewalk_test,return
elseif strcmpi(varargin{1}, '--f')
    type makefigs_ridgewalk
    makefigs_ridgewalk;
    return
end

if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else
    dt=1;
end

%Default values of parameters
alpha=1/4;
alg='amp';
N=1.5;
chi=0;
fmax=[];
fmin=[];
params=[];
derivparams=[];

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='col';
end

if iscell(varargin{end})
    params=varargin{end};
    varargin=varargin(1:end-1);
end

if allall(size(varargin{1})==size(varargin{end}))
    mask=varargin{end};
    varargin=varargin(1:end-1);    
else
    mask=[];
end

fs=varargin{end};
varargin=varargin(1:end-1);

%/********************************************************************
%Sorting out input params
if ~isempty(params)
    if ischar(params{end})
        alg=params{end};
        params=params(1:end-1);
    end
    if length(params)>=4
        fmax=params{1}(:).*dt;
        fmin=params{2}(:).*dt;
        params=params(3:end);
    end
    if length(params)>=1
        N=params{1};
    end
    if length(params)>=2
        chi=params{2};
    end
    if length(params)>=3
        alpha=params{3};
    end
end
%\********************************************************************
%fmin,fmax,N,chi

%/********************************************************************
if length(fs)<3
    %Contingency for short input frequency array
    for i=1:nargout
        varargout{i}=[];
    end
    return
    disp('Sorry, RIDGEWALK needs at least three frequencies in the wavelet')
    disp('transform in order to identify ridges.')
end
%\********************************************************************


if length(varargin)>1
    w=zeros([size(varargin{1},1) size(varargin{1},2) length(varargin)]);
    for i=1:length(varargin);
        w(:,:,i,:)=varargin{i};
    end
else
    w=varargin{1};
end

if size(w,3)
    disp(['RIDGEWALK detecting a set of ' int2str(size(w,3)) ' transforms.'])
end

[bool,rq,wjoint,omjoint]=isridgepoint(w,fs,chi,alg,fmin,fmax,mask);

%figure,jpcolor(omjoint'),
%figure,jpcolor(diff(omjoint)')

disp('RIDGEWALK chaining ridges.')
[id,ir,jr,wr,fr]=ridgechains(fs,N,bool,wjoint,omjoint,alpha,mask);
%figure,plot(diff(fr))

for i=1:nargout
    varargout{i}=[];
end

%/*************************************************************************
%Compute bias parameters
br=[];
cr=[];
if ~isempty(id)
    [id,ir,jr]=colbreaks(id,ir,jr);
    
    if strcmpi(str,'mat')
        na=nargout-4+1;
    elseif strcmpi(str,'col')
        na=(nargout-2)./size(w,3)-2+1;
    end
    if na<=0
        [a,om]=instmom(dt,w,'endpoint');
        [wr,fr]=ridgeinterp(fs,rq,ir,jr,w,om);
    elseif na>0
        %Output bandwidth & curvature or deviation vectors
        if size(w,3)==1
            %Just normal bandwith and curvature if it's a single time series
            [a,om,up,curv]=instmom(dt,w,'endpoint');
            [wr,fr,br,cr]=ridgeinterp(fs,rq,ir,jr,w,om,up,curv);
        else
            %Deviation vectors if it's a set of time series
            w1=vdiff(w,1,'endpoint')./dt;
            w2=vdiff(w1,1,'nans')./dt;
            [a,om]=instmom(dt,w,'endpoint');
            [wr,fr,w1r,w2r]=ridgeinterp(fs,rq,ir,jr,w,om,w1,w2);
            fbar=vrep(vmean(fr,2,squared(wr)),size(wr,2),2);
            wrmag=vrep(sqrt(sum(squared(wr),2)),size(wr,2),2);
            br=(w1r-sqrt(-1)*fbar.*wr)./wrmag;                  %See (17) of L012
            cr=(w2r-2*sqrt(-1)*fbar.*w1r-fbar.^2.*wr)./wrmag;   %See (18) of L012
            %Fix for extra nans due to derivative
            cr(1,:)=cr(2,:);
        end
    end
end
%\*************************************************************************


%/*************************************************************************
%Re-organizing output

%Using the fact that COL2CELL works for matrices, and for empties
%vsize(ir,jr,wr,fr,br,cr)
col2cell(ir,jr,wr,fr,br,cr);

varargout{1}=ir;
varargout{2}=jr;
if ~isempty(id)
    if strcmpi(str,'mat')
        varargout{3}=wr;
        varargout{4}=fr;
        if na>=1
            varargout{5}=br;
        end
        if na>=2
            varargout{6}=cr;
        end
    elseif strcmpi(str,'col')
        %Put into one-column cell array format
        N=size(w,3);
        varargout(3:3+N-1)=splitcells(wr);
        varargout(3+N:3+2*N-1)=splitcells(fr);
        if ~isempty(br)&&na>=1
            varargout(3+2*N:3+3*N-1)=splitcells(br);
        end
        if ~isempty(cr)&&na>=2
            varargout(3+3*N:3+4*N-1)=splitcells(cr);
        end
    end
end
%\*************************************************************************


disp('RIDGEWALK finished.')
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

ga=2;
be=4;

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,ga,be,fs,'bandpass'},'mirror');
[wp,wn]=vectmult(tmat,wx,wy);

%Form ridges of component time series
[ir,jr,wr,fr]=ridgewalk(dt,wn,fs,{1.5,0,'phase'}); 
[ir2,jr2,wr2,fr2]=ridgewalk(dt,wn,fs,{1.5,0,'amplitude'});

err=vsum(abs(wr{1}-wr2{1}).^2,1)./vsum(abs(wr{1}).^2,1);
reporttest('RIDGEWALK phase and amplitude signal estimation error for NPG-06',err<1e-3)

[ir,jr,wpr,wnr]=ridgewalk(dt,wp,wn,fs,{3,0});   
[ir2,jr2,wxr,wyr]=ridgewalk(dt,wx,wy,fs,{3,0});   
[wpr2,wnr2]=vectmult(tmat,wxr{1},wyr{1});

err=vmean((abs([wpr{1} wnr{1}]-[wpr2 wnr2])./sqrt(abs([wpr{1} wnr{1}]).^2+abs([wpr2 wnr2])).^2),1);
reporttest('RIDGEWALK XY vs. PN invariance joint ridge for NPG-06',allall(err<1e-10))

%[ir,jr,wr,fr]=ridgewalk(dt,wn(:,1),fs(1),true(size(wn(:,1))),{1.5,0,'phase'}); 

% [ir,jr,wr,fr,br,cr]=ridgewalk(dt,wp,wn,fs,{3,0},'mat');   
% 
% %Check to see if these mean what I think they mean
% wr=wr{1};
% fr=fr{1};
% br=br{1};
% cr=cr{1};
% 
% om=frac(sum(squared(wr).*fr,2),sum(squared(wr),2));
% om(:,2)=om;
% 
% br2=vdiff(vdiff(wr,1),1)./dt./dt-2*1i*vdiff(wr,1).*om./dt-om.^2.*wr;
% br2=br2./abs(wr);
% 
% 




%[ir,jr,wr]=ridgewalk(dt,wp,wn,fs,{3,0},'mat');   
%[ir,jr,wpr,wnr]=ridgewalk(dt,wp,wn,fs,{3,0});   
