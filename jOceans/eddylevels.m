function[rho,ymid,S,Snoise,Srat]=eddylevels(varargin)
%EDDYLEVELS  Eddy ridge significance levels using the survival function.
%
%   RHO=EDDYLEVELS(DX,XMID,YBIN,X,Y,XNOISE,YNOISE,N) returns the eddy ridge
%   significance levels RHO for a Lagrangian dataset, based on comparision 
%   of the survival function of the data with that of an ensemble of noise. 
%
%   For details and theory on the use of EDDYLEVELS, see
%
%       Lilly, J. M. and P. Perez-Brunius (2021). Extracting statistically
%           significant eddy signals from large Lagrangian datasets using
%           wavelet ridge analysis, with application to the Gulf of Mexico.
%           Submitted to Nonlinear Processes in Geohysics.
%
%   This paper will be referred to as LPB21 afterwards.
%
%   The eddy ridges for the data and the noise are assumed to have been 
%   calculated using EDDYRIDGES, while the noise itself is generated using 
%   NOISEDRIFTERS.  See MAKE_GOMED for an example.
%
%   The significance levels are computed in terms of X and Y, two ridge-
%   averaged quantities.  XNOISE and YNOISE are the comparable quantities 
%   in the noise ensemble, which is N times the size of the data. 
%
%   In LPB21, the ridge-averaged nondimensional frequency OMEGA_AST_BAR is
%   used for X, and the ridge length L times the square of the circularity 
%   XI is used for Y.  See therein for a discussion and other possibilties.
%
%   A survival function is created for Y as a function of X, with the 
%   x-bins centered on midpoints XMID and spanning XMID-DX/2 to XMID+DX/2,
%   and with YBIN giving the bin edges for the Y variable.  
%
%   RHO is then the significance level of the ridges based on interpolating
%   within the ratio of the two survival functions.  Smaller values of RHO 
%   correspond to higher significance levels.  RHO has the size of X and Y. 
%
%   [RHO,YMID,S,SNOISE,SRAT]=EDDYLEVELS(...) also returns the y-bin 
%   midpoints YMID, the data and noise survival functions S and SNOISE, 
%   and their ratio SRAT=SNOISE./S.  The survival functions are all 
%   matrices with LENGTH(YMID) rows and LENGTH(XMID) columns.  
%
%   EDDYLEVELS(...,'symmetric') optionally forces the noise survival
%   function to be symmetric about its x-axis, as expected if the x-axis
%   spans a symmetric band of frequencies centered about zero as in LPB21.
%
%   EDDYLEVELS(...,'sort') optionally sorts the survival function ratio
%   SRAT along its y-axis in descending order, in order to remove
%   small-scale jaggedness that can arise in sparsely sampled bins.  This
%   is the recommended setting, for reasons discussed in LPB21.
%
%   Usage: rho=eddylevels(dx,xmid,ybin,x,y,xnoise,ynoise,N);
%          [rho,ymid,S,Snoise,Srat]
%                 =eddylevels(dx,xmid,ybin,x,y,xnoise,ynoise,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2021 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    eddylevels_test,return
end
 
dx=varargin{1};
xmid=varargin{2};
ybins=varargin{3};
x=varargin{4};
y=varargin{5};
xnoise=varargin{6};
ynoise=varargin{7};
N=varargin{8};

symstr='nosym';
sortstr='nosort';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'sym')||strcmpi(varargin{end}(1:3),'nosym')
            symstr=varargin{end};
        elseif  strcmpi(varargin{end}(1:3),'sor')||strcmpi(varargin{end}(1:3),'nosor')
            sortstr=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

mat=zeros(length(ybins)-1,length(xmid));
for i=1:length(xmid)
      bool=(x>xmid(i)-dx/2)&(x<xmid(i)+dx/2);
      [mat(:,i),~,ymid]=twodhist(x(bool),y(bool),[-inf inf],ybins);
end
S=flipud(cumsum(flipud(mat),1));

matnoise=zeros(length(ybins)-1,length(xmid));
for i=1:length(xmid)
      bool=(xnoise>xmid(i)-dx/2)&(xnoise<xmid(i)+dx/2);
      matnoise(:,i)=twodhist(xnoise(bool),ynoise(bool),[-inf inf],ybins);
end
Snoise=flipud(cumsum(flipud(matnoise),1))./N;

if strcmpi(symstr(1:3),'sym')
   %force this to be exactly symmetric, as there is no reason it should not be
   Snoise=frac(1,2)*(Snoise+fliplr(Snoise));
end

Srat=Snoise./S;

if strcmpi(sortstr(1:3),'sor')
    Srat_sorted=Srat;
    Srat_sorted(isinf(Srat))=0;
    Srat_sorted(isnan(Srat))=0;
    Srat=sort(Srat_sorted,'descend');
end

rho=10.^interp2(xmid,ymid,log10(Srat),x,y);


