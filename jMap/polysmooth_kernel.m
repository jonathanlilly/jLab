function[W]=polysmooth_kernel(dist,ws,B,kernstr)
%POLYSMOOTH_KERNEL  Returns the weighting kernel employed by POLYSMOOTH.
%
%   POLYSMOOTH_KERNEL is a low-level function called by POLYSMOOTH. 
%
%   W=POLYSMOOTH_KERNEL(DS,[],B,KERNSTR) where DS is an array of distances 
%   and B is a bandwidth, returns the value of the weighting function 
%   specified by the string KERNSTR at distances DS.  
%
%   See POLYSMOOTH for valid choices of the string KERNSTR.
%
%   B is either a scalar, or a matrix the same size as DS(:,:,1).
%
%   W=POLYSMOOTH_KERNEL(DS,[],B,L) where L is a scalar uses a Gaussian with 
%   a standard deviation set to B/L, as in W=EXP(-1/2*(L*DS/B)^2).
%
%   W=POLYSMOOTH_KERNEL(DS,[],B,K) where K is an array takes K to be a 
%   custom kernel, with K(1) corresponding to DS/B=0 and K(end) to DS/B=1, 
%   interpolates within K to find the weight W at any distance.
%
%   W=POLYSMOOTH_KERNEL(T,[],TAU,KERNSTR) works in just the same way if T 
%   is a time deviation and TAU is a temporal bandwidth.
%
%   W=POLYSMOOTH_KERNEL(DS,WS,B,KERNSTR) optionally multiplies W, after 
%   computed as described above, by the additional weights given in WS.  
%
%   See also POLYSMOOTH.
%
%   Usage: W=polysmooth_kernel(d,[],B,kernstr);
%          W=polysmooth_kernel(d,ws,B,kernstr);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018--2020 J.M. Lilly --- type 'help jlab_license' for details

if length(B(:))~=1
    B=vrep(B,size(dist,3),3);
end

dist=frac(dist,B);
W=zeros(size(dist));
bool=(dist<=1);
dist=dist(bool);

if ischar(kernstr)
    if strcmpi(kernstr(1:3),'epa')
        W(bool)=frac(2,pi).*(1-dist.^2);
    elseif strcmpi(kernstr(1:3),'gau')
        W(bool)=frac(9,2*pi*(1-exp(-9/2))).*exp(-frac((3*dist).^2,2));
    elseif strcmpi(kernstr(1:3),'bis')
        W(bool)=frac(3,pi).*(1-dist.^2).^2;
    elseif strcmpi(kernstr(1:3),'tri')
        W(bool)=frac(220,81*pi).*(1-dist.^3).^3;
    elseif strcmpi(kernstr(1:3),'uni')
        W(bool)=frac(1,pi);
    end
elseif length(kernstr)==1
    %Variable Gaussian kernel
%    W(bool)=exp(-(1/2)*(dist.*kernstr).^2)-exp(-(1/2)*(kernstr).^2);
    W(bool)=frac(kernstr.^2,2*pi*(1-exp(-kernstr.^2/2))).*exp(-(1/2)*(dist.*kernstr).^2);
else
    %Custom kernel
    K=kernstr(:);
    W(bool)=interp1((0:length(K)-1)'./(length(K)-1),K,dist);
end

if ~isempty(ws)
    W=W.*ws;
end

%elseif iscell(kernstr)
%   if strcmpi(kernstr{1}(1:3),'hyb')
%   %Hybrid Gaussian kernel
%          W(bool)=exp(-(1/2)*(dist.*kernstr).^2)-exp(-(1/2)*(kernstr).^2);
%   elseif strcmpi(kernstr{1}(1:3),'mor')
%   %Morse kernel
%g=exp(-(1/2)*(x./(1/3)).^2);
% x=[0:0.01:1];
% sigma=1./[1:.5:10]';
% for i=1:length(sigma)
%     G(:,i)=exp(-(1/2)*(x./sigma(i)).^2)-exp(-(1/2)*(1./sigma(i)).^2);
%     G(:,i)=G(:,i)./max(G(:,i));
% end
%figure,plot(x,G)

% if strcmpi(kernstr(1:3),'epa')
%     %W=frac(2,pi).*(1-dist.^2).*(1+sign(1-dist));
%     W=(1-dist.^2).*(1+sign(1-dist));
% elseif strcmpi(kernstr(1:3),'gau')
%     %W=frac(2,pi.*B.^2).*exp(-frac(dist.^2,2)).*(1+sign(1-dist./N));
%     W=exp(-frac((3*dist).^2,2)).*(1+sign(1-dist));
% elseif strcmpi(kernstr(1:3),'bis')
%     W=((1-dist.^2).^2).*(1+sign(1-dist));
% elseif strcmpi(kernstr(1:3),'tri')
%     W=((1-dist.^3).^3).*(1+sign(1-dist));
% end

W=vswap(W,nan,0);
