function[varargout]=ridgeinterp(varargin)
% RIDGEINTERP  Interpolate quantity values onto ridge locations.
%
%   RIDGEINTERP is a low-level function called by RIDGEWALK.
%
%   XI=RIDGEINTERP(FS,RQ,IR,JR,X) where RQ is a "ridge quantity" at
%   *radian* frequencies FS, and IR and JR are time and scale indices of 
%   ridges, interpolates quantity X along the ridge locations to give XI.
%
%   IR and JR give indices into the first two dimensions of X.  The output
%   XI has the the same number of rows and columns as IR and JR, that is,
%   SIZE(XI,1)=SIZE(IR,1) and SIZE(XI,2)=SIZE(IR,2).  The locations of 
%   NANs in IR and JR are duplicated in XI. 
%
%   X may have more than one 'page' along its third dimension, in which 
%   case each page is interpolated separately, and SIZE(XI,3)=SIZE(X,3). 
%
%   RQ is output by RIDGEQUANTITY based on a wavelet transform output by
%   WAVETRANS, and IR and JR are output by RIDGEWALK.  
%
%   RIDGEINTERP interpolates transform values between discrete frequency 
%   levels to find a more precise value of the transform along the ridges
%   than simply looking up the values of X at rows IR and columns JR.  
%
%   XI=RIDGEINTERP(FS,RQ,IR,JR,MU,X) interprets the ridges as belonging
%   to the MU-th derivative of the signal to which the quantity X belongs
%   and applies an appropriate correction factor. 
%
%   [XI1,XI2,...,XIN]=RIDGEINTERP(FS,RQ,IR,JR,X1,X2,...,XN) also 
%   interpolates the quantities X1,X2,...,XN, all having the same number of 
%   rows and columns as IR and JR.
%
%   RIDGEINTERP uses fast quadratic interpolation via QUADINTERP.  In rare
%   cases, quadratic interpolation fails for an individual point and 
%   therefore linear interpolation is used instead.
%   __________________________________________________________________
%
%   RIDGEINTERP is a low-level function called by RIDGEWALK.
%
%   See also RIDGEWALK, RIDGEQUANTITY, QUADINTERP.
%
%   Usage:  xi=ridgeinterp(rq,ir,jr,x);
%           [xi1,xi2,xi3]=ridgeinterp(rq,ir,jr,x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2018 J.M. Lilly --- type 'help jlab_license' for details    


rq=varargin{1};
ir=varargin{2};
jr=varargin{3};
varargin=varargin(4:end);

if iscell(varargin{end})||isempty(varargin{end})
    derivparams=varargin{end};
    varargin=varargin(1:end-1);
else
    derivparams=[];
end

index=find(~isnan(ir));
if ~isempty(index)
    varargout=ridgeinterp1_quadratic(index,ir,jr,derivparams,rq,varargin);
else
    for i=1:length(varargin)
        varargout{i}=[];
    end
end
              
function[outargs]=ridgeinterp1_quadratic(index,ir,jr,derivparams,rq,args)
sizeir=size(ir);
vcolon(ir,jr);
vindex(ir,jr,index,1);


indexr=nonnan(sub2ind(size(rq),ir,jr));

%Ridge quantity along the ridges, and at one scale up and down

di=size(rq,1);
dr=rq(indexr);
drp=rq(indexr+di);
drn=rq(indexr-di);

[~,jre]=quadinterp(jr-1,jr,jr+1,abs(drn).^2,abs(dr).^2,abs(drp).^2);

%Rare complete failure of quadratic interpolation associated 
%with the ridge quantity changing sign between jrp and jrn, yet
%not having the ridge quantity at jr being between these two values
        
bool=~( (jr+1>jre) & (jre> jr-1)); 
jre(bool)=lininterp(jr(bool)-1,jr(bool)+1,drn(bool),drp(bool));

if ~isempty(derivparams)
    ga=derivparams{1};
    be=derivparams{2};
    mu=derivparams{3};
    
    fact=frac(morsefreq(ga,be-mu),morsefreq(ga,be));
    
    jro=frac(jr,fact);
    jrp=frac(jr+1,fact);
    jrn=frac(jr-1,fact);
    
    btop=(ceil(jro)>=(size(rq,2)-1)); 
    jro(btop)=size(rq,2)-1;
    jrp(btop)=size(rq,2);
    jrn(btop)=size(rq,2)-2;
    
    bbottom=(floor(jro)<=2); 
    jro(bbottom)=2;
    jrp(bbottom)=3;
    jrn(bbottom)=1;

    %The problem here is that I also have to interpolate to find the
    %values of the quantities to be interpolated at the rescaled scale
    %levels, so you have to do another linear interpolation between 
    %the floors and the ceilings.  You can't just do the above any
    %more because the rescaled scale levels are no longer jr+/-1.
    
%     figure,plot(jro)
%     maxmax(jro)
%     minmin(jro)
    indexrf=nonnan(sub2ind(size(rq),ir,floor(jro)));
    indexrc=nonnan(sub2ind(size(rq),ir,ceil(jro)));
    
    indexrpf=nonnan(sub2ind(size(rq),ir,floor(jrp)));
    indexrpc=nonnan(sub2ind(size(rq),ir,ceil(jrp)));
    
    indexrnf=nonnan(sub2ind(size(rq),ir,floor(jrn)));
    indexrnc=nonnan(sub2ind(size(rq),ir,ceil(jrn)));
end


for i=1:length(args)
   x=args{i};
   if isreal(x)
        xr=nan*ones(sizeir(1),size(x,3));
   else
        xr=(nan+sqrt(-1)*nan)*ones(sizeir(1),size(x,3));
   end
   for k=1:size(x,3)
        xk=x(:,:,k);
        if isempty(derivparams)
            xro=xk(indexr);
            xrp=xk(indexr+di);
            xrn=xk(indexr-di);      
        else
            xro=lininterp(floor(jro),ceil(jro),xk(indexrf),xk(indexrc),jro); 
            xrp=lininterp(floor(jrp),ceil(jrp),xk(indexrpf),xk(indexrpc),jrp); 
            xrn=lininterp(floor(jrn),ceil(jrn),xk(indexrnf),xk(indexrnc),jrn); 
        end

        xrk=quadinterp(jr-1,jr,jr+1,xrn,xro,xrp,jre);         
        %Use linear interpolation where quadratic fails
        xrk(bool)=lininterp(jr(bool)-1,jr(bool)+1,xrn(bool),xrp(bool),jre(bool)); 
        xr(index,k)=xrk;
   end
   %Linearly interpolate between the approximate ridge and the bracketing curve 
   outargs{i}=xr;
end

