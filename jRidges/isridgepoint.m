function[bool,rq,w,om]=isridgepoint(w,fs,delta,str,fmin,fmax,mask)
%ISRIDGEPOINT  Finds wavelet ridge points using one of several criterion.
%
%   ISRIDGEPOINT is a low-level function called by RIDGEWALK.
%  
%   BOOL=ISRIDGEPOINT(W,FS,DELTA,STR) where W is a wavelet transform matrix  
%   at *radian* frequecies FS, finds all ridge points of W with amplitudes
%   |W| exceeding the amplitude cutoff DELTA.  Several different different 
%   ridge defintions may be used and are specified by STR.
%
%   BOOL is a matrix of the same size as W, which is equal to one for 
%   those elements of W which are ridge points, and zero otherwise.
%
%   STR may be either of the following:
%
%        'phase'       Rate of transform change of phase definition
%        'amplitude'   Maxima of transfom amplitude definition
%
%   For all definitions, ISRIDGEPOINT rejects spurious ridge points.
%   These tend to occur on the flanks of interesting signals, and 
%   reflect the wavelet structure rather than the signal structure.
%
%   A ridge point is considered spurious if either it is located at an
%   amplitude minima, or if the frequency anomaly (transform frequency
%   minus scale frequency) is a maximum.
%
%   BOOL=ISRIDGEPOINT(W,FS,DELTA,STR,FMIN,FMAX) only returns ridge points
%   between frequencies FMIN and FMAX, which may be either scalars or 
%   arrays of length SIZE(W,1).
%
%   BOOL=ISRIDGEPOINT(W,FS,DELTA,STR,FMIN,FMAX,MASK), where MASK is a boolean
%   variable of the same size as W, only returns ridge points at locations
%   at which MASK is true.
%
%   [BOOL,RQ,W,OMEGA]=ISRIDGEPOINT(...) also returns the ridge quantity RQ,
%   the complex-valued joint wavelet transform W, and the transform's
%   instantaneous frequency OMEGA, all of which are of the same size. 
%
%   See also RIDGEINTERP, RIDGEWALK.
%
%   Usage: [bool,rq,w,om]=isridgepoint(w,fs,delta,str);
%          [bool,rq,w,om]=isridgepoint(w,fs,delta,str,fmin,fmax);
%          [bool,rq,w,om]=isridgepoint(w,fs,delta,str,fmin,fmax,mask);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2020 J.M. Lilly --- type 'help jlab_license' for details
 
%        'groove'      Joint amplitude / phase definition



disp('RIDGEWALK looking for ridge points...')

%Re-doing way of handling amplitude ridges to remove ridge breaking
if size(w,3)~=1
    [a,om]=instmom(w,1,3);
else
    [a,om]=instmom(w);
end

if strcmpi(str(1:3),'amp')
    rq=a;
elseif strcmpi(str(1:3),'pha')
    fsmat=vrep(fs(:)',size(w,1),1);
    rq=om-fsmat;
end    
if size(w,3)~=1
    phaseavg=frac(sum(abs(w).*w,3),sum(abs(w).^2,3));
    w=sqrt(sum(abs(w).^2,3)).*rot(angle(phaseavg));
end

%figure,jpcolor(rq'),shading flat

rqm=circshift(rq,+1,2);
rqp=circshift(rq,-1,2);

if strcmpi(str(1:3),'amp')
   bool=(rqm<=rq)&(rqp<=rq);
elseif strcmpi(str(1:3),'pha')
   %This is d/ds < 0 since scale decreases in columns
   %actually it's d/dln s
   bool=(rqm<0&rqp>=0)|(rqm<=0&rqp>0);
end
%[ii,jj]=find(bool);
%figure,plot(ii,jj,'.')

%d/dlns=s d/ds since  dlns = 1/s ds 
%Does d/dlns =0 where d/ds=0?  Seems so

err=abs(rq);

%Ensure maximum not minimum
bool((bool&circshift(bool,-1,2))&err>circshift(err,-1,2))=0; 
bool((bool&circshift(bool,+1,2))&err>circshift(err,+1,2))=0; 

bool1= ~isnan(w);   %Remove NANs
bool2=~(abs(w)<delta);  %Remove those less than cutoff amplitude
bool=bool.*bool1.*bool2;
bool(:,[1 end])=0;

%Running FMIN and FMAX 
if ~isempty(fmin)
    %figure,plot(fsmat)
    %figure,plot(om)
    %    figure,plot(om)fmin,fmax
    %    hlines(fmin),hlines(fmax)
        
    fmin=vrep(fmin,size(w,2),2);
    fmax=vrep(fmax,size(w,2),2);
    if size(fmin,1)==1
        fmin=vrep(fmin,size(w,1),1);
        fmax=vrep(fmax,size(w,1),1);
    end
    %vsize(fmin,fmax,fsmat)
    
    %This way is more intuitive, and seems better at rejecting short
    %ridges.  Otherwise the results can be identical.
    bool3=(om>fmin)&(om<fmax);
    %bool3=(fsmat>=fmin)&(fsmat<=fmax)&(om>fmin)&(om<fmax);
    %fsmat=vrep(fs(:)',size(w,1),1);
    %bool3=(fsmat>=fmin)&(fsmat<=fmax);
    bool=bool.*bool3;
    %figure,plot(om),hlines(fs)
end

if ~isempty(mask)
    bool=bool.*mask;
end
disp(['RIDGEWALK found ' int2str(length(find(bool))) ' ridge points.'])
