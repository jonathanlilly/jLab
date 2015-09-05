function[varargout]=morsederiv(nmax,ga,be,om)
%MORSEDERIV  Frequency-domain derivatives of generalized Morse wavelets.
%
%   DCELL=MORSEDERIV(N,GAMMA,BETA,OMEGA) returns the first N normalized
%   frequency derivatives of the lowest-order generalized Morse wavelet 
%   specified by parameters GAMMA and BETA, at the radian frequency OMEGA.
%
%   DCELL is a cell array with DCELL{N} containing the Nth derivative.
%
%   This derivative is defined as
%
%           dn = ( omega^n / psi ) * ( d^n psi / d omega^n )
%
%   where omega is the radian frequency and psi is the wavelet.
% 
%   N must be a scalar. The other input parameters must either be matrices
%   of the same size, or some may be matrices and the others scalars.
%
%   DCELL=MORSEDERIV(N,GAMMA,BETA) evaluates the derivative at the peak
%   frequency by default.
%
%   For details see
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   See also MORSEWAVE, MORSEMOM, BELLPOLY.
%
%   'morsewave --t' runs a test.
%
%   Usage: dcell=morsederiv(4,gamma,beta,omega);
%          dcell=morsederiv(4,gamma,beta);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2009 J.M. Lilly --- type 'help jlab_license' for details  

if strcmpi(nmax,'--t')
  morsederiv_test;return
end

if nargin==3
    om=morsefreq(ga,be);
end

% %Earlier algorithm using bellpoly
% for n=1:nmax;
%     kcell{n}=morsekappa(ga,be,om,n);
%     psi{n}=om.^n.*bellpoly(kcell(1:n));
% end

%Quicker to use cum2mom directly
kcell{1}=0*om;
for n=1:nmax;
    kcell{n+1}=morsekappa(ga,be,om,n);
end
psi=cum2mom(kcell);
psi=psi(2:end);


if nargout==1&&nmax>1
   varargout{1}=psi;
else
   varargout=psi;
end    
    
function[kappan]=morsekappa(ga,be,om,n)
%Derivatives of natural log of wavelet
%See Lilly and Olhede (2009)

prodn=ones(size(ga));
for p=0:n-1
    prodn=prodn.*(ga-p);
end

kappan=be.*(-1).^(n-1).*factorial(n-1)-om.^ga.*prodn;


function[psi1,psi2,psi3,psi4]=morsederiv_peak(ga,be)
%Algebraic expressions at peak frequency
%See Lilly and Olhede (2009)

[psi1,psi2,psi3,psi4]=vzeros(size(ga));
psi2=-be.*ga;
psi3=-be.*ga.*(ga-3);
psi4=(3.*(be.*ga).^2-be.*(6+(ga-1).*(ga-2).*(ga-3)));


function[]=morsederiv_test

ga1=(1/2:.1:10);
be1=(1/2:.1:10);
[ga,be]=meshgrid(ga1,be1);
om=morsefreq(ga,be);

[psi1a,psi2a,psi3a,psi4a]=morsederiv(4,ga,be,om);
[psi1b,psi2b,psi3b,psi4b]=morsederiv_peak(ga,be);

tol=1e-6;
bool(1)=aresame(psi1a,psi1b,tol);
bool(2)=aresame(psi2a,psi2b,tol);
bool(3)=aresame(psi3a,psi3b,tol);
bool(4)=aresame(psi4a,psi4b,tol);
reporttest('MORSEDERIV peak expressions',all(bool))


