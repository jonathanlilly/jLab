function[varargout]=windtrans(varargin)
%WINDTRANS  Ekman-like transfer-functions for the wind-driven response.
%
%   G=WINDTRANS(OMEGA,Z,FC,DELTA,MU,H) returns the no-slip transfer 
%   function for the wind-driven currents evaluated at frequencies OMEGA
%   and at depths Z.  G will have LENGTH(OMEGA) rows and LENGTH(Z) columns.
%
%   Here FC is the local Coriolis frequency, DELTA is the Ekman depth.
%   MU is the Madsen depth, and H is the boundary layer depth.
%
%   The units of these quantitites are important.  For consistency with 
%   other routines, OMEGA and FC are in radians per day, while Z, DELTA,
%   MU, and H are all in meters.  The units of G are then m^2 s / kg.
%
%   WINDTRANS(...,'free') instead returns the free slip transfer function.
%
%   For details on the expressions implemented by this function, see
%
%        Lilly, J. M. and S. Elipot (2021). A unifying perspective on 
%            transfer function solutions to the unsteady Ekman problem. 
%            Fluids, 6 (2): 85, 1--36. 
%   __________________________________________________________________
%
%   Special forms
%
%   By default WINDTRANS uses the general transfer function.  WINDTRANS 
%   will also employ limiting expressions for special cases, as follows.
%
%       H = Inf                --  Mixed Ekman / Madsen solution
%       H = Inf, MU = 0        --  Ekman solution
%       H = Inf, DELTA = 0     --  Madsen solution
%       MU = 0                 --  Finite-layer Ekman solution
%       DELTA = 0              --  Finite-layer Madsen solution
%
%   These have to be coded separately because the full solution is singular
%   in these cases.
%   __________________________________________________________________
%
%   Computational options
%
%   When computed over a wide range of parameter space, the transfer
%   function tends to be encounter numerical overflow when the arguments to
%   Bessel functions become large, causing its computation to fail.
%
%   To avoid this problems, by default WINDTRANS switches to using a highly
%   accurate thirty-term expansion about the large-argument exponential 
%   behavior of the Bessel functions when their arguments exceed 10^2.9. 
%
%   Two other options are available, primarily for testing purposes. Both
%   of these other algorithms lead to artifacts and are not recommended.
%
%   WINDTRANS(...'far',...) switches instead to use the (inferior) one-term 
%   expansion, also known as the far-inertial limit.
%
%   WINDTRANS(...,'general',...) uses the general formula with no switch.
%
%   Note that these options only apply to no-slip solution.  The free-slip
%   solution, which is not deemed to be physically relevant, is only
%   computed with the general formula.
%
%   For details on these algorithms, see Lilly and Elipot (2021).
%   __________________________________________________________________
%
%   'windtrans --t' runs a some tests.
%
%   Usage: G=windtrans(omega,z,fc,delta,mu,h);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019--2021 J.M. Lilly --- type 'help jlab_license' for details
 

%   WINDTRANS(OMEGA,Z,FC,DELTA,MU,A,'AMP') alternately parameterizes the
%   transfer function in terms of the magntide of the near-inertial peak A
%   in 1/m, together with the Ekman and Madsen depths DELTA and MU.

%   Gradient output
%
%   [G,DDELTA,DH]=WINDTRANS(OMEGA,Z,FC,DELTA,0,H) for the Ekman-like model
%   also returns the partial derivatives with respect to DELTA and H.

if strcmp(varargin{1}, '--t')
    windtrans_test,return
end
 
omega=varargin{1}(:)/24/3600;%convert to rad/second
z=varargin{2};
fc=varargin{3}/24/3600;%convert to rad/second
delta=varargin{4};
mu=varargin{5};
h=varargin{6};

tol=100;
rho=1027;
str='expansion';%use the 30-term tilde expansion by default
%str='two';%use the two-term tilde expansion by default
slipstr='noslip';

bool=false(2,1);
for i=7:nargin
    istr=varargin{i}(1:3);
    if strcmpi(istr,'gen')||strcmpi(istr,'eli')||strcmpi(istr,'far')||strcmpi(istr,'ada')||strcmpi(istr,'exp')||strcmpi(istr,'two')
        str=varargin{i};
        if strcmpi(str(1:3),'exp')
            bool(1)=true;
        end
    elseif strcmpi(istr,'fre')||strcmpi(istr,'nos')
        slipstr=varargin{i};
        if strcmpi(slipstr(1:3),'fre')
            bool(2)=true;
        end
    end
end

%this is just to keep track of an unsuitable combination of input arguments
if all(bool)
     disp('Sorry, the tilde expansion algorithm is not implemented for the free-slip transfer function.')  
end
%these may have been either input, or reflecting the default behavior
if strcmpi(str(1:3),'exp')&&strcmpi(slipstr(1:3),'fre')
     str='gen';
end

omega1=omega;
if length(z)~=1
    [z,omega]=meshgrid(z,omega);
end

%if isinf(h)||delta==0||mu==0
%    str='lilly';%use the lilly version of expressions for special cases
%end

%str,slipstr
G=zeros(size(omega));
if strcmpi(slipstr(1:3),'nos')
    %----------------------------------------------------------------------
    %The noslip forms
        %[G,xiz,xih,xi0]=windtrans_general_noslip(delta,fc,rho,z,mu,omega,h,str);
   %     G=windtrans_general_noslip(delta,fc,rho,z,mu,omega,h,str);
    if strcmpi(str(1:3),'gen')|| strcmpi(str(1:3),'far')||strcmpi(str(1:3),'exp')||strcmpi(str(1:3),'two')
        [G,varargout{2},varargout{3}]=windtrans_lilly_noslip(delta,fc,rho,z,mu,omega,h,str);
    elseif strcmpi(str(1:3),'eli')%asymptotic forms, elipot edition
        G=windtrans_elipot_noslip(delta,fc,rho,z,mu,omega,h);
    end
elseif strcmpi(slipstr(1:3),'fre')
    %----------------------------------------------------------------------
    %The freeslip forms

    K0=frac(1,2)*delta.^2.*abs(fc);
    K1=frac(1,2)*mu*abs(fc);

    s=sign(fc).*sign(1+omega./fc);
    zo=delta.^2./mu;
    [xiz,xih,xi0]=xis(s,zo,delta,z,omega,fc,h);

    coeff=frac(sqrt(2).*rot(-s.*pi/4),delta.*rho.*abs(fc).*sqrt(abs(1+omega./fc)));
    [k0z,i0z,k1h,i1h,k10,i10]=bessels_freeslip(xiz,xih,xi0);
 
    if strcmpi(str(1:3),'gen')||(strcmpi(str(1:3),'lil')&&((K0~=0)&&(K1~=0)))
        numer=i0z.*k1h+i1h.*k0z;
        denom=i1h.*k10-i10.*k1h;
        G=coeff.*frac(numer,denom);
    elseif strcmpi(str(1:3),'lil')
        %    [h,K0,K1]
        if mu==0
            coeff=frac(sqrt(2).*rot(-s*pi/4),delta.*abs(fc).*rho.*sqrt(abs((1+omega./fc))));
            cosharg=sqrt(2).*rot(s*pi/4).*frac(h-z,delta).*sqrt(abs((1+omega./fc)));
            sinharg=sqrt(2).*rot(s*pi/4).*frac(h,delta).*sqrt(abs((1+omega./fc)));
            G=coeff.*frac(cosh(cosharg),sinh(sinharg));
        elseif delta==0
            coeff=frac(2,rho*K1);
            k0z=besselk(0,2*rot(s*pi/4).*sqrt(frac(z,K1./abs(fc)).*abs(1+omega./fc)));
            k1h=besselk(1,2*rot(s*pi/4).*sqrt(frac(h,K1./abs(fc)).*abs(1+omega./fc)));
            i0z=besseli(0,2*rot(s*pi/4).*sqrt(frac(z,K1./abs(fc)).*abs(1+omega./fc)));
            i1h=besseli(1,2*rot(s*pi/4).*sqrt(frac(h,K1./abs(fc)).*abs(1+omega./fc)));
            G=coeff.*(k0z+frac(k1h.*i0z,i1h));
        end
    elseif strcmpi(str(1:3),'eli')
        
        delta1=sqrt(2.*K0./(omega+fc));
        delta2=K1./(omega+fc);
        xize=2*sqrt(1i*(zo+z)./delta2);
        xi0e=2*sqrt(1i*zo./delta2);
        xihe=2*sqrt(1i*(zo+h)./delta2);
        
        if mu==0
            coeff=frac(1,rho.*sqrt(1i*(omega+fc).*K0));
            numer=cosh((1+1i).*(h-z)./delta1);
            denom=sinh((1+1i).*h./delta1);
            G=coeff.*numer./denom;
        elseif delta==0
            coeff=(2./rho./K1);
            k0z=besselk(0,2.*sqrt(1i.*z./delta2));
            k1h=besselk(1,2*sqrt(1i.*h./delta2));
            i0z=besseli(0,2*sqrt(1i.*z./delta2));
            i1h=besseli(1,2*sqrt(1i.*h./delta2));
            G=coeff.*(k0z+frac(k1h.*i0z,i1h));
        else
            coeff=frac(1,rho.*sqrt(1i.*(omega+fc).*K0));
            k0z=besselk(0,xize);
            k10=besselk(1,xi0e);
            i10=besseli(1,xi0e);
            i0z=besseli(0,xize);
            k1h=besselk(1,xihe);
            i1h=besseli(1,xihe);
            numer=i0z.*k1h+i1h.*k0z;
            denom=i1h.*k10-i10.*k1h;
            G=coeff.*frac(numer,denom);
        end
    end
end
G(z>h)=nan;

varargout{1}=G;
%varargout{2}=xih;
%varargout{3}=xiz;
%varargout{4}=xi0;

function[xiz,xih,xi0]=xis(s,zo,delta,z,omega,fc,h)

xiz=2*sqrt(2).*rot(s.*pi/4).*(zo./delta).*sqrt((1+z./zo).*abs((1+omega./fc)));
xih=2*sqrt(2).*rot(s.*pi/4).*(zo./delta).*sqrt((1+h./zo).*abs((1+omega./fc)));
xi0=2*sqrt(2).*rot(s.*pi/4).*(zo./delta).*sqrt(abs((1+omega./fc)));

function[k0z,i0z,k0h,i0h,k10,i10]=bessels_noslip(varargin)

argz=varargin{1};
argh=varargin{2};
if length(varargin)==3
    arg0=varargin{3};
end

k0z=besselk(0,argz);
i0z=besseli(0,argz);
k0h=besselk(0,argh);
i0h=besseli(0,argh);

if nargout>4
    k10=besselk(1,arg0);
    i10=besseli(1,arg0);
end

function[k0z,i0z,k0h,i0h,k10,i10]=besseltildes_noslip(argz,argh,arg0,nterms)

k0z=besselktilde(0,argz,nterms);
i0z=besselitilde(0,argz,nterms);
k0h=besselktilde(0,argh,nterms);
i0h=besselitilde(0,argh,nterms);

k10=besselktilde(1,arg0,nterms);
i10=besselitilde(1,arg0,nterms);

function[k0z,i0z,k1h,i1h,k10,i10]=bessels_freeslip(argz,argh,arg0)

k0z=besselk(0,argz);
i0z=besseli(0,argz);
k1h=besselk(1,argh);
i1h=besseli(1,argh);

if nargout>4
    k10=besselk(1,arg0);
    i10=besseli(1,arg0);
end

function[G]=windtrans_expansion_noslip(delta,fc,rho,z,mu,omega,h,str)

zo=delta.^2./mu;
s=sign(fc).*sign(1+omega./fc);
[xiz,xih,xi0]=xis(s,zo,delta,z,omega,fc,h);
coeff=frac(sqrt(2).*rot(-s.*pi/4),delta.*rho.*abs(fc).*sqrt(abs(1+omega./fc)));

if strcmpi(str(1:3),'two')
    [k0z,i0z,k0h,i0h,k10,i10]=besseltildes_noslip(xiz,xih,xi0,2);
elseif strcmpi(str(1:3),'exp')
    [k0z,i0z,k0h,i0h,k10,i10]=besseltildes_noslip(xiz,xih,xi0,30);
end
        
numer=exp(xi0-xiz).*i0h.*k0z-exp(xi0+xiz-2*xih).*k0h.*i0z;
denom=i0h.*k10+exp(2*xi0-2*xih).*k0h.*i10;
G=coeff.*frac(numer,denom);

bool=(omega==-fc);
if ~isinf(h)&&~isinf(zo)
    if length(z)==1
        G(bool)=frac(4*zo,rho*abs(fc)*delta.^2).*frac(sqrt(1+h./zo)-sqrt(1+z./zo),(1+z./zo).^(1/4));
    else
        G(bool)=frac(4*zo,rho*abs(fc)*delta.^2).*frac(sqrt(1+h./zo)-sqrt(1+z(bool)./zo),(1+z(bool)./zo).^(1/4));
    end
elseif isinf(h)
    G(bool)=inf;
end

function[G,xiz,xih,xi0]=windtrans_general_noslip(delta,fc,rho,z,mu,omega,h,str)

zo=delta.^2./mu;
s=sign(fc).*sign(1+omega./fc);
[xiz,xih,xi0]=xis(s,zo,delta,z,omega,fc,h);
coeff=frac(sqrt(2).*rot(-s.*pi/4),delta.*rho.*abs(fc).*sqrt(abs(1+omega./fc)));

[k0z,i0z,k0h,i0h,k10,i10]=bessels_noslip(xiz,xih,xi0);

%numerically, better to do it like this to avoid having huge numbers dominate
%tests shows this removes the large-argument overflow when only h (not z) is large

numer=k0z./k10-(k0h./k10).*(i0z./i0h);
denom=1+(k0h./k10).*(i10./i0h);
G=coeff.*frac(numer,denom);


if strcmpi(str(1:3),'exp')||strcmpi(str(1:3),'two')
    bool=log10(abs(xiz))>2.9;
    if ~isempty(find(bool,1)) 
        if length(z)==1
            G(bool)=windtrans_expansion_noslip(delta,fc,rho,z,mu,omega(bool),h,str);
        else
            G(bool)=windtrans_expansion_noslip(delta,fc,rho,z(bool),mu,omega(bool),h,str);
        end
    end
elseif strcmpi(str(1:3),'far')
    bool=log10(abs(xiz))>2.9;
    if ~isempty(find(bool,1)) 
        if length(z)==1
            G(bool)=windtrans_farinertial_noslip(delta,fc,rho,z,mu,omega(bool),h);
        else
            G(bool)=windtrans_farinertial_noslip(delta,fc,rho,z(bool),mu,omega(bool),h);
        end
    end
end
G=windtrans_inertiallimit(G,delta,fc,rho,z,mu,omega,h);




function[G]=windtrans_farinertial_noslip(delta,fc,rho,z,mu,omega,h)

zo=delta.^2./mu;
s=sign(fc).*sign(1+omega./fc);
[xiz,xih,xi0]=xis(s,zo,delta,z,omega,fc,h);
coeff=frac(sqrt(2).*rot(-s.*pi/4),delta.*abs(fc).*rho);

numer=exp(-xiz+xi0)-exp(xiz+xi0-2*xih);
denom=1+exp(2*xi0-2*xih);
G=coeff.*frac(1,sqrt(abs(1+omega./fc)).*(1+z./zo).^(1/4)).*frac(numer,denom);
%Don't apply the inertial limit!

function[G]=windtrans_inertiallimit(G,delta,fc,rho,z,mu,omega,h)

zo=delta.^2./mu;
bool=(omega==-fc);

if ~isempty(find(bool,1))
    if ~isinf(h)&&~isinf(zo)
        if length(z)==1
            %G(bool)=frac(2*zo,rho*abs(fc)*delta.^2).*log(frac(1+h./zo,1+z./zo));
            G(bool)=frac(2,rho*abs(fc)*mu).*log(frac(1+h./zo,1+z./zo));
        else
            %G(bool)=frac(2*zo,rho*abs(fc)*delta.^2).*log(frac(1+h./zo,1+z./zo));
            G(bool)=frac(2,rho*abs(fc)*mu).*log(frac(1+h./zo,1+z(bool)./zo));
        end
    elseif ~isinf(h)&&isinf(zo)
        if length(z)==1
            G(bool)=frac(2,rho*abs(fc)*delta.^2).*(h-z);
        else
            G(bool)=frac(2,rho*abs(fc)*delta.^2).*(h-z(bool));
        end
    else
        G(bool)=inf;
    end
end

function[G,ddelta,dh]=windtrans_lilly_noslip(delta,fc,rho,z,mu,omega,h,str)

zo=delta.^2./mu;
s=sign(fc).*sign(1+omega./fc);
[xiz,xih,xi0]=xis(s,zo,delta,z,omega,fc,h);
%[delta,zo,K0,K1]

if h==inf
    if mu==0
        %Ekman solution
        coeff=frac(sqrt(2).*rot(-s.*pi/4),delta.*abs(fc).*rho);
        G=coeff.*frac(exp(-(1+s.*1i).*(z./delta).*sqrt(abs((1+omega./fc)))),sqrt(abs((1+omega./fc))));
    elseif delta==0
        %Madsen solution
        coeff=frac(4,rho*abs(fc)*mu);
        G=coeff.*besselk(0,2*sqrt(2)*rot(s.*pi/4).*sqrt(frac(z,mu).*abs(1+omega./fc)));
    else
        %Mixed solution
        k0z=besselk(0,xiz);
        k10=besselk(1,xi0);
        coeff=frac(sqrt(2).*rot(-s.*pi/4),delta.*abs(fc).*rho.*sqrt(abs((1+omega./fc))));
        G=coeff.*frac(k0z,k10);
    end
else
    if mu==0
        %finite-layer Ekman
        coeff=frac(sqrt(2).*rot(-s*pi/4),delta.*abs(fc).*rho.*sqrt(abs((1+omega./fc))));
        %sinharg=sqrt(2).*rot(s*pi/4).*frac(h-z,delta).*sqrt(abs((1+omega./fc)));
        %cosharg=sqrt(2).*rot(s*pi/4).*frac(h,delta).*sqrt(abs((1+omega./fc)));
        %G=coeff.*frac(sinh(sinharg),cosh(cosharg));    
        %can't do it this way because you divide two huge numbers
        
        argh=sqrt(2).*rot(s*pi/4).*frac(h,delta).*sqrt(abs((1+omega./fc)));
        argz=sqrt(2).*rot(s*pi/4).*frac(z,delta).*sqrt(abs((1+omega./fc)));
        
        numer=exp(-argz)-exp(argz).*exp(-2*argh);
        denom=1+exp(-2*argh);
        G=coeff.*frac(numer,denom);
        
        bool=(omega==-fc);
        if length(z)==1
            G(bool)=frac(2,rho*abs(fc)*delta.^2).*(h-z);
        else
            G(bool)=frac(2,rho*abs(fc)*delta.^2).*(h-z(bool));
        end
    elseif delta==0
        %         if z==0
        %             coeff=frac(4,rho*abs(fc)*mu);
        %             argz=2*sqrt(2)*rot(s*pi/4).*sqrt(frac(z,mu).*abs(1+omega./fc));
        %             argh=2*sqrt(2)*rot(s*pi/4).*sqrt(frac(h,mu).*abs(1+omega./fc));
        %
        %             [k0z,i0z,k0h,i0h]=bessels_noslip(argz,argh);
        %             length(find(isnan(k0z)))
        %             length(find(isnan(i0z)))
        %             length(find(isnan(k0h)))
        %             length(find(isnan(i0h)))
        %             length(find(isnan(argh)))
        %
        %             %figure,plot(abs(argh))
        %
        %             G=-coeff.*(log(xi0)+frac(1,2)*(z./zo)+frac(k0h,i0h));
        %             G=-coeff.*(log(xi0./xih)+frac(1,2)*(z./zo));
        %             length(find(isnan(xi0./xih)))
        %
        %         else
        %finite-layer Madsen
        coeff=frac(4,rho*abs(fc)*mu);
        argz=2*sqrt(2)*rot(s*pi/4).*sqrt(frac(z,mu).*abs(1+omega./fc));
        argh=2*sqrt(2)*rot(s*pi/4).*sqrt(frac(h,mu).*abs(1+omega./fc));
        [k0z,i0z,k0h,i0h]=bessels_noslip(argz,argh);
        G=coeff.*(k0z-i0z.*frac(k0h,i0h));
        
        bool=(omega==-fc);
        if length(z)==1
            G(bool)=frac(1,2)*coeff.*log(h./z);
        else
            G(bool)=frac(1,2)*coeff.*log(h./z(bool));
        end
        %end
    else
        %General solution
        G=windtrans_general_noslip(delta,fc,rho,z,mu,omega,h,str);
    end
end

if mu==0
    s=sign(fc).*sign(1+omega./fc);
    Gamma=sqrt(2).*rot(s*pi/4).*sqrt(abs((1+omega./fc)));
    ddelta1=(Gamma.*frac(h,delta).*tanh(Gamma.*frac(h,delta))-1).*G./delta;
    numer=exp(Gamma.*frac(-z,delta))+exp(-Gamma.*frac(2*h-z,delta));
    denom=1+exp(-Gamma.*frac(2*h,delta));
    ddelta2=-frac(2,delta.^2.*abs(fc).*rho).*frac(h-z,delta).*frac(numer,denom);
    ddelta=ddelta1+ddelta2;
    dh1=-Gamma.*frac(1,delta).*tanh(Gamma.*frac(h,delta)).*G;
    dh2=frac(2,delta.^2.*abs(fc).*rho).*frac(numer,denom);
    dh=dh1+dh2;
else
    ddelta=[];
    dh=[];
end


function[G]=windtrans_elipot_noslip(delta,fc,rho,z,mu,omega,h)

zo=delta.^2./mu;
K0=frac(1,2)*delta.^2.*abs(fc);
K1=frac(1,2)*mu*abs(fc);

delta1=sqrt(2.*K0./(omega+fc));
delta2=K1./(omega+fc);
xiz=2*sqrt(1i*(zo+z)./delta2);
xi0=2*sqrt(1i*zo./delta2);
xih=2*sqrt(1i*(zo+h)./delta2);

if h==inf
    if K1==0
       %Ekman solution
        coeff=frac(1,rho.*sqrt(1i*(omega+fc).*K0));
        G=coeff.*exp(-z.*(1+1i)./delta1);
    elseif K0==0
        %Madsen solution
        coeff=(2./rho./K1);
        k0z=besselk(0,2.*sqrt(1i.*z./delta2));
        G=coeff.*k0z;
    else
        %Mixed solution
        coeff=frac(1,rho.*sqrt(1i.*(omega+fc).*K0));
        G=coeff.*besselk(0,xiz)./besselk(1,xi0);
    end
else
    if K1==0
        %finite-layer Ekman
        coeff=frac(1,rho.*sqrt(1i*(omega+fc).*K0));
        numer=sinh((1+1i).*(h-z)./delta1);
        denom=cosh((1+1i).*h./delta1);
        G=coeff.*numer./denom;
    elseif K0==0
        %finite-layer Madsen
        coeff=(2./rho./K1);
        k0z=besselk(0,2.*sqrt(1i.*z./delta2));
        k0h=besselk(0,2*sqrt(1i.*h./delta2));
        i0z=besseli(0,2*sqrt(1i.*z./delta2));
        i0h=besseli(0,2*sqrt(1i.*h./delta2));
        G=coeff.*(k0z-frac(k0h.*i0z,i0h));
    else
        %General solution
        coeff=frac(1,rho.*sqrt(1i.*(omega+fc).*K0));
        [k0z,i0z,k0h,i0h,k10,i10]=bessels_noslip(xiz,xih,xi0);
        numer=i0h.*k0z-k0h.*i0z;
        denom=i10.*k0h+k10.*i0h;
        G=coeff.*frac(numer,denom);
    end
end

function[]=windtrans_test
windtrans_test_gradient;
windtrans_test_limits;

function[]=windtrans_test_gradient

delta=10.^(-1:0.05:3);
h=10.^(log10(15.15):.05:5);
[dg,hg]=meshgrid(delta,h);

[G,G1,G2,G3,G4,dg1,dg2]=vzeros(size(dg));
ddelta=1e-6;dh=1e-6;
omega=-1e-4;
%omega=-1.5e-4;
%omega=-5e-4;
for j=1:length(delta)
    for i=1:length(h)
        %G=windtrans(omega,z,fc,delta,mu,h);
       [G(i,j),dg1(i,j),dg2(i,j)]=windtrans(omega,15,1e-4,dg(i,j),0,hg(i,j));
        G1(i,j)=windtrans(omega,15,1e-4,dg(i,j)-ddelta/2,0,hg(i,j));
        G2(i,j)=windtrans(omega,15,1e-4,dg(i,j)+ddelta/2,0,hg(i,j));
        G3(i,j)=windtrans(omega,15,1e-4,dg(i,j),0,hg(i,j)-dh/2);
        G4(i,j)=windtrans(omega,15,1e-4,dg(i,j),0,hg(i,j)+dh/2);
    end
end

dg1hat=frac(G2-G1,ddelta);
dg2hat=frac(G4-G3,dh);

bool=dg1hat~=0;
eps1=maxmax(log10(frac(abs(dg1hat(bool)-dg1(bool)),sqrt(squared(dg1hat(bool))+squared(dg1(bool))))));
bool=dg2hat~=0;
eps2=maxmax(log10(frac(abs(dg2hat(bool)-dg2(bool)),sqrt(squared(dg2hat(bool))+squared(dg2(bool))))));

%jpcolor(log10(delta),log10(h),real(dg1))
%figure,jpcolor(log10(delta),log10(h),log10(frac(abs(dg1hat-dg1),sqrt(squared(dg1hat)+squared(dg1)))))
%jpcolor(log10(delta),log10(h),log10(frac(abs(dg2hat-dg2),sqrt(squared(dg2hat)+squared(dg2)))))

reporttest(['WINDTRANS analytic and numerical gradients match for Ekman case'],eps1<-4&&eps2<-4)

function[]=windtrans_test_limits
 
N=1000;
z=[0.1:1:100];
h=[inf 200];
K0=[0 1/10 1/10];
K1=[1 0 1];
fc=1e-4;
delta=sqrt(frac(2*K0,fc));
mu=frac(2*K1,fc);

[bool1,bool2,bool3,bool4]=vzeros(length(K0),2);
[Ge,Go,Gl,Ga]=vzeros(N,length(z),length(K0),2);
slipstr='noslip';
omega=fourier(1,N,'two');
           
for s=[1 -1]
    for i=1:length(K0)
        for j=1:2
            Ge(:,:,i,j)=windtrans(omega,z,s/2,delta(i),mu(i),h(j),'elipot',slipstr);
            Gl(:,:,i,j)=windtrans(omega,z,s/2,delta(i),mu(i),h(j),slipstr);
            Go(:,:,i,j)=windtrans(omega,z,s/2,delta(i)+1e-14,mu(i)+0.5,h(j),slipstr,'general');
            Ga(:,:,i,j)=windtrans(omega,z,s/2,delta(i)+1e-14,mu(i)+0.5,h(j),slipstr,'two');
            Ge(:,:,i,j)=Ge(:,:,i,j)./maxmax(abs(Gl(:,:,i,j)));
            Go(:,:,i,j)=Go(:,:,i,j)./maxmax(abs(Gl(:,:,i,j)));
            Ga(:,:,i,j)=Ga(:,:,i,j)./maxmax(abs(Gl(:,:,i,j)));
            Gl(:,:,i,j)=Gl(:,:,i,j)./maxmax(abs(Gl(:,:,i,j)));
            
            bool1(i,j)=aresame(Ge(:,:,i,j),Gl(:,:,i,j),1e-8);
            bool2(i,j)=aresame(Ge(:,:,i,j),Go(:,:,i,j),1e-1);
            bool3(i,j)=aresame(Gl(:,:,i,j),Go(:,:,i,j),1e-1);
            bool4(i,j)=aresame(Ga(:,:,i,j),Go(:,:,i,j),1e-20);
        end
    end
    
    if s==1
        reporttest(['WINDTRANS Elipot and Lilly forms match for f>0 and ' slipstr ' condition'],allall(bool1))
        reporttest(['WINDTRANS general and special forms match for f>0 and ' slipstr ' condition'],allall(bool2).*allall(bool3))
        reporttest(['WINDTRANS general and expansion forms match for f>0 and ' slipstr ' condition'],allall(bool4).*allall(bool3))
    else
        reporttest(['WINDTRANS Elipot and Lilly forms match for f<0 and ' slipstr ' condition'],allall(bool1))
        reporttest(['WINDTRANS general and special forms match for f<0 and ' slipstr ' condition'],allall(bool2).*allall(bool3))
        reporttest(['WINDTRANS general and expansion forms match for f>0 and ' slipstr ' condition'],allall(bool4).*allall(bool3))
    end
end


% older test for freeslip forms; not currently  particularly relevant 
% h=20;
% z=[0.1:1:20];
% [bool1,bool2,bool3]=vzeros(length(K0),1);
% [Ge,Go,Gl]=vzeros(N,length(z),length(K0));
% slipstr='freeslip';
% for s=[1 -1]
%     for i=1:length(K0)
%         for j=1:length(h)
%             i,j,s
%             Ge(:,:,i,j)=windtrans(omega,z,s/2,delta(i),mu(i),h(j),'elipot',slipstr);
%             Gl(:,:,i,j)=windtrans(omega,z,s/2,delta(i),mu(i),h(j),'lilly',slipstr);
%             Go(:,:,i,j)=windtrans(omega,z,s/2,delta(i)+1e-10,mu(i)+0.5,h(j),slipstr,'general');
%             Ge(:,:,i,j)=Ge(:,:,i,j)./maxmax(abs(Gl(:,:,i,j)));
%             Go(:,:,i,j)=Go(:,:,i,j)./maxmax(abs(Gl(:,:,i,j)));
%             Gl(:,:,i,j)=Gl(:,:,i,j)./maxmax(abs(Gl(:,:,i,j)));
%             bool1(i,j)=aresame(Ge(:,:,i,j),Gl(:,:,i,j),1e-8);
%             bool2(i,j)=aresame(Ge(:,:,i,j),Go(:,:,i,j),1e-1);
%             bool3(i,j)=aresame(Gl(:,:,i,j),Go(:,:,i,j),1e-1);
%         end
%     end
%     
%     if s==1
%         reporttest(['WINDTRANS Elipot and Lilly forms match for f>0 and ' slipstr ' condition'],allall(bool1))
%         reporttest(['WINDTRANS general and special forms match for f>0 and ' slipstr ' condition'],allall(bool2).*allall(bool3))
%     else
%         reporttest(['WINDTRANS Elipot and Lilly forms match for f<0 and ' slipstr ' condition'],allall(bool1))
%         reporttest(['WINDTRANS general and special forms match for f<0 and ' slipstr ' condition'],allall(bool2).*allall(bool3))
%     end
% end


