function[W]=equivkern(varargin)
%EQUIVKERN  Equivalent kernels for local polynomial fitting.
%
%   WS=EQUIVKERN(DS,XS,YS,RHO,H,ALPHA,BETA) computes the equivalent kernel
%   WS for a order RHO local polynomial fit using the generalized beta 
%   kernel with parameters ALPHA and BETA, and bandwidth H.
% 
%   DS, XS, and YS are cell arrays containing the distance, x deviations, 
%   and y deviations between observation points and mapping points as 
%   output by those functions.  WS is a cell arry of the same size as DS.
%
%   EQUIVKERN is to be called after calling TWODSORT or SPHERESORT and 
%   before calling POLYSMOOTH.  
% 
%   For RHO=1 or RHO=2, using the equivalent kernel can greatly accelarate 
%   POLYFIT for problems in which the same design (the distribution of 
%   observation points) is used multiple times. 
%
%   To use with POLYSMOOTH, use WS as a weight function in the input to 
%   POLYSMOOTH, and call POLYSMOOTH with a uniform kernel (BETA=0) and 
%   bandwidth H. See POLYSMOOTH for details.
%
%   The idea of equivalent kernels goes back to Silverman (1984).  The 
%   expressions used herein for linear and quadratic fits are taken from
% 
%       Takeda, Farsiu, and Milanfar (2007), "Kernel regression for image 
%           processing and reconstruction", IEEE Transactions on Image 
%           Processing, v. 16 (2), 349--366.
%
%   _________________________________________________________________
%
%   XXX dumped from POLYMAP 
%   Acceleration for repeated designs
%
%   For problems in which observations are carried out multiple times using
%   the same design, or distribution of observational points, computations
%   of the field value (but not, at the moment, its derivatives) for RHO>0 
%   can be greatly accerated via the use of equivalent kernels.
% 
%   The idea behind the equivalent kernel is that a fit to a polynomial
%   results in a weighted sum of the data points.  This implies that a 
%   higher-order fit be expressed as a zeroth-order fit using a data-
%   depedent kernel called the equivalent kernel.
% 
%   The equivalent kernel depends on the original choice of smoothing
%   kernel and bandwidth, together with the observational design.
%   Importantly, it does not depend upon the values of the observations.
%
%   Thus, we can compute the equivalent kernel one time for, say, a
%   first-order fit.  Then we can generate many maps with different
%   observed data values but the same design, implementing the first-order 
%   fit through its equivalent, but faster, zeroth-order representation.
%
%   To make use of the equivalent kernel, after calling TWODSORT or 
%   SPHERESORT, one would then call EQUIVKERN followed by POLYMAP,
%
%        WS=EQUIVKERN(DS,XS,YS,RHO,H,ALPHA,BETA);
%        ZHAT=POLYMAP(DS,XS,YS,[],ZS,WS,RHO,{H,ALPHA,0});
%  
%   where in the input to POLYMAP, the weights WS output by EQUIVKERN
%   have the effect of weighting the data points.  Importantly, in the call
%   to POLYMAP, one should set BETA to zero to specify a uniform kernel.
%   This has the effect of not modifying the equivalent kernel weights WS.
%
%   After the equivalent kernel is computed, POLYMAP proceeds about 4x
%   faster for a linear fit, or 10x faster for a quadratic fit.
%
%   The idea of equivalent kernels goes back to Silverman (1984).  The 
%   expressions used herein for linear and quadratic fits are taken from
% 
%       Takeda, Farsiu, and Milanfar (2007), "Kernel regression for image 
%           processing and reconstruction", IEEE Transactions on Image 
%           Processing, v. 16 (2), 349--366.
%
%   While it is possible to formulate the equivalent kernels for the
%   derivatives of the filed, in BETA, this has not yet been implemented.
%
%   See POLYSMOOTH for more details.
%
%   'equivkern --t' runs a test.
%
%   Usage: W=equivkern(ds,xs,ys,rho,H,alpha,beta);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    equivkern_test,return
end
 
%DS,XS,YS
ds=varargin{1};
xs=varargin{2};
ys=varargin{3};
rho=varargin{4};
H=varargin{5};
al=varargin{6};
be=varargin{7};

if length(al)==1
    al=al+zeros(length(ds),size(ds{1},2));
end
if length(be)==1
    be=be+zeros(length(ds),size(ds{1},2));
end
if length(H)==1
    H=H+zeros(length(ds),size(ds{1},2));
end

if rho==0
    disp('You don''t need an equivalent kernel for a constant fit.')
    W=[];
    return
end

W=equivkern_direct_linear(ds,xs,ys,rho,H,al,be);

function[W]=equivkern_direct_linear(ds,xs,ys,rho,H,al,be)

kappa=frac(1,2*pi)*frac(al,beta(2./al,be+1));

for i=1:length(ds)
    %plvsize(ds{i},H(i,:))
    Wi=kappa(i,:).*(1-(ds{i}./H(i,:)).^al(i,:)).^be(i,:);

    xvect=zeros(2,1,size(ds{i},1),size(ds{i},2));
    xvect(1,1,:,:)=xs{i};
    xvect(2,1,:,:)=ys{i};

    s0(1,1,:,:)=vsum(Wi,1);

    svect=zeros(2,1,1,size(ds{i},2));
    S=zeros(2,2,1,size(ds{i},2));

    svect(1,:,:)=vsum(xs{i}.*Wi,1);
    svect(2,:,:)=vsum(ys{i}.*Wi,1);
    S(1,1,:,:)=vsum(xs{i}.*xs{i}.*Wi,1);
    S(1,2,:,:)=vsum(xs{i}.*ys{i}.*Wi,1);
    S(2,1,:,:)=S(1,2,:,:);
    S(2,2,:,:)=vsum(ys{i}.*ys{i}.*Wi,1);

    Sinv=matinv(S);

    if rho==1
        svectTSinv=pagemtimes(svect,'transpose',Sinv,'none');
        numer=1-pagemtimes(svectTSinv,xvect);
        denom=s0-pagemtimes(svectTSinv,svect);
        %vsize(s0,S,xvect,svect,numer,denom,Wi)
        W{i,1}=squeeze(frac(numer,denom)).*Wi;
        %W{i,1}=Wi;
    elseif rho==2
        s13=zeros(1,3,1,size(ds{i},2));
        s23=zeros(2,3,1,size(ds{i},2));
        s33=zeros(3,3,1,size(ds{i},2));

        s13(1,1,:,:)=S(1,1,:,:);
        s13(1,2,:,:)=S(2,1,:,:);
        s13(1,3,:,:)=S(2,2,:,:);
 
        %this matrix has only four unique entries
        s23(1,1,:,:)=vsum(xs{i}.*xs{i}.*xs{i}.*Wi,1);
        s23(1,2,:,:)=vsum(xs{i}.*xs{i}.*ys{i}.*Wi,1);
        s23(1,3,:,:)=vsum(xs{i}.*ys{i}.*ys{i}.*Wi,1);
        s23(2,1,:,:)=s23(1,2,:,:);
        %s23(2,1,:,:)=vsum(ys{i}.*xs{i}.*xs{i}.*Wi,1);
        s23(2,2,:,:)=s23(1,3,:,:);
        %s23(2,2,:,:)=vsum(ys{i}.*xs{i}.*ys{i}.*Wi,1);
        s23(2,3,:,:)=vsum(ys{i}.*ys{i}.*ys{i}.*Wi,1);

        %this matrix has only five unique entries
        s33(1,1,:,:)=vsum(xs{i}.*xs{i}.*xs{i}.*xs{i}.*Wi,1);
        s33(1,2,:,:)=vsum(xs{i}.*xs{i}.*xs{i}.*ys{i}.*Wi,1);
        s33(1,3,:,:)=vsum(xs{i}.*xs{i}.*ys{i}.*ys{i}.*Wi,1);
        s33(2,1,:,:)=s33(1,2,:,:);%from symmetry
        %s33(2,1,:,:)=vsum(xs{i}.*ys{i}.*xs{i}.*xs{i}.*Wi,1);
        s33(2,2,:,:)=s33(1,3,:,:);
        %s33(2,2,:,:)=vsum(xs{i}.*ys{i}.*xs{i}.*ys{i}.*Wi,1);
        s33(2,3,:,:)=vsum(xs{i}.*ys{i}.*ys{i}.*ys{i}.*Wi,1);
        s33(3,1,:,:)=s33(1,3,:,:);%from symmetry
        s33(3,2,:,:)=s33(2,3,:,:);%from symmetry
        %s33(3,1,:,:)=vsum(ys{i}.*ys{i}.*xs{i}.*xs{i}.*Wi,1);
        %s33(3,2,:,:)=vsum(ys{i}.*ys{i}.*xs{i}.*ys{i}.*Wi,1);
        s33(3,3,:,:)=vsum(ys{i}.*ys{i}.*ys{i}.*ys{i}.*Wi,1);

        s33inv=matinv(s33);
        
        %s22=S
        %svect'=s12

        S12=permute(svect,[2 1 3 4])-pagemtimes(pagemtimes(s13,s33inv),'none',s23,'transpose');
        S22=S-pagemtimes(pagemtimes(s23,s33inv),'none',s23,'transpose');
        S13=s13-pagemtimes(pagemtimes(svect,'transpose',Sinv,'none'),s23);
        S33=s33-pagemtimes(pagemtimes(s23,'transpose',Sinv,'none'),s23);

        S22inv=matinv(S22);
        S33inv=matinv(S33);

        xvech=zeros(3,1,size(ds{i},1),size(ds{i},2));
        xvech(1,1,:,:)=xs{i}.*xs{i};
        xvech(2,1,:,:)=xs{i}.*ys{i};
        xvech(3,1,:,:)=ys{i}.*ys{i};

        S12S22inv=pagemtimes(S12,S22inv);
        S13S33inv=pagemtimes(S13,S33inv);

        numera=pagemtimes(S12S22inv,xvect);
        numerb=pagemtimes(S13S33inv,xvech);
        denoma=pagemtimes(S12S22inv,svect);
        denomb=pagemtimes(S13S33inv,'none',s13,'transpose');
        W{i,1}=squeeze(frac(1-numera-numerb,s0-denoma-denomb)).*Wi;
    end
end



function[W]=equivkern_loop(ds,xs,ys,H,al,be)

kappa=frac(1,2*pi)*frac(al,beta(2./al,be+1));

% for i=1:length(ds)
%     for j=1:size(ds{i},2)
% 
%         Wij=kappa(i,j).*(1-(ds{i}(:,j)./H(i,j)).^al(i,j)).^be(i,j);
% 
%         
%         xvect=zeros(2,1,size(ds{i},1),size(ds{i},2));
%         s0(1,1,:,:)=vsum(Wi,1);
% 
%         svect=zeros(2,1,1,size(ds{i},2));
%         S=zeros(2,2,1,size(ds{i},2));
% 
%         svect(1,:,:)=vsum(xs{i}.*Wi,1);
%         svect(2,:,:)=vsum(ys{i}.*Wi,1);
%         S(1,1,:,:)=vsum(xs{i}.*xs{i}.*Wi,1);
%         S(1,2,:,:)=vsum(xs{i}.*ys{i}.*Wi,1);
%         S(2,1,:,:)=S(1,2,:,:);
%         S(2,2,:,:)=vsum(ys{i}.*ys{i}.*Wi,1);
% 
%         Sinv=matinv(S);
%         svectTSinv=pagemtimes(svect,'transpose',Sinv,'none');
%         numer=1-pagemtimes(svectTSinv,xvect);
%         denom=s0-pagemtimes(svectTSinv,svect);
%         %vsize(s0,S,xvect,svect,numer,denom,Wi)
%         W{i,1}=squeeze(frac(numer,denom)).*Wi;
%         %W{i,1}=Wi;
%     end
% end



function[]=equivkern_test
 
[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:600);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.125:3);
yo=(-3:.125:3);

B=1;

[ds,xs,ys,zs]=twodsort(xdata,ydata,zdata,xo,yo,B);
tic;z1=polysmooth(ds,xs,ys,[],zs,[],1,{B,2,1});etime1=toc;
ws=equivkern(ds,xs,ys,1,B,2,1);
tic;z1b=polysmooth(ds,xs,ys,[],zs,ws,0,{B,2,0});etime2=toc;

reporttest('EQUIVKERN first order fit ',maxmax(abs(z1-z1b)<1e-10))
disp(['POLYSMOOTH was ' num2str(etime1./etime2) ' times faster using EQUIVKERN.'])
[ds,xs,ys,zs]=twodsort(xdata,ydata,zdata,xo,yo,B);

tic;z2=polysmooth(ds,xs,ys,[],zs,[],2,{B,2,1});etime1=toc;
ws=equivkern(ds,xs,ys,2,B,2,1);
tic;z2b=polysmooth(ds,xs,ys,[],zs,ws,0,{B,2,0});etime2=toc;

reporttest('EQUIVKERN second order fit',maxmax(abs(z2-z2b)<1e-10))
disp(['POLYSMOOTH was ' num2str(etime1./etime2) ' times faster using EQUIVKERN.'])


