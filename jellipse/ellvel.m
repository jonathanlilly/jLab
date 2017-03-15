function[v]=ellvel(varargin)
% ELLVEL  Average and instantaneous ellipse velocities. 
%
%   V=ELLVEL(KAPPA,LAMBDA,THETA,PHI,STR) where KAPPA and LAMBDA are the 
%   amplitude and linearity of a time-varying ellise, THETA is its time-
%   varying orientation, and PHI is its time-varying phase, returns various
%   measures of the ellispe 'velocity'.
%
%   STR determines the velocity quantity to be output.
%
%       STR       Symbol     Description
%       ----------------------------------------------------------------
%      'geo'      VM         Geometric mean velocity    
%      'ave'      VBAR       Period-averaged speed          
%      'kin'      VEKE       Kinetic energy velocity        
%      'cir'      VGAMMA     Circulation velocity        
%      'azi'      VAZ        Instantaneous azimuthal velocity        
%      'ins'      VI         Instantaneous speed 
%
%   All are signed quantities reflecting the direction of motion.  That is, 
%   a velocity is defined to be positive when the ellipse is orbited in the 
%   mathematically positive (counterclockwise) sense, and negative when the
%   ellipse is orbited in the mathematically negative sense.
%
%   VM=ELLVEL(KAPPA,LAMBDA,THETA,PHI) with STR omitted returns VM, the
%   geometric mean velocity.  This is a basic of ellipse velocity that is 
%   analagous to the geometric mean radius as a measure of ellipse size.
%
%   ELLVEL(DT,...,) optionally uses DT as the data sample rate, with a 
%   default value of DT=1.  DT is a scalar.
%
%   ELLVEL(DT,...,FACT,STR) optionally converts the physical units of 
%   velocity through a multiplication by FACT, with a default value of 
%   FACT=1.  For example, FACT=1e5 converts kilometers into centimeters, 
%   while FACT=1e5/24/3600 converts km/day into cm/sec. 
%   ____________________________________________________________________
%   
%   Cell array input/output
%
%   If ELLVEL is given cell array input, it returns cell array output.
%
%   Thus KAPPA, LAMBDA, THETA, and PHI may each be cell arrays of the 
%   same size, where each element in the cell array is a numerical array.
%
%   VM and other velocity measures will be also cell arrays of this size.
%
%   In this case ELLVEL(DT,...) also works with DT a scalar or an 
%   array whose length is the number of elements in the cell arrays.
%   ____________________________________________________________________
%
%   See also ELLRAD, ELLPARAMS, ELLDIFF.
%
%   Usage:  v=ellvel(kappa,lambda,theta,phi);
%           v=ellvel(dt,kappa,lambda,theta,phi,1e5);
%           vbar=ellvel(kappa,lambda,theta,phi,'ave');
%
%   'ellvel --f' generates a sample figure.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2014 J.M. Lilly --- type 'help jlab_license' for details    


if strcmpi(varargin{1}, '--f')
  type makefigs_ellvel 
  makefigs_ellvel
  return
end

dt=1;
if iscell(varargin{2})
    if ~iscell(varargin{1})
        dt=varargin{1};
        varargin=varargin(2:end);
    end
elseif length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
end

if ischar(varargin{end})
    str=varargin{end}(1:3);
    varargin=varargin(1:end-1);
else
    str='geo';
end

if length(varargin{end})==1
    fact=varargin{end};
    varargin=varargin(1:end-1);
else 
    fact=1;
end

kappa=varargin{1};
lambda=varargin{2};
theta=varargin{3};
phi=varargin{4};

if ~isempty(kappa)
    if ~iscell(kappa)
        v=ellvel_one(dt,kappa,lambda,theta,phi,fact,str);
    else
        v=kappa;
        for i=1:length(kappa)
            if isscalar(dt)
                dti=dt;
            else
                dti=dt(i);
            end
            v{i}=ellvel_one(dti,kappa{i},lambda{i},theta{i},phi{i},fact,str);
        end
    end
else
    v=kappa;
end


function[v]=ellvel_one(dt,kappa,lambda,theta,phi,fact,str)

strend='endpoint';
vswap(lambda,0,1e-10);
vswap(lambda,1,1-1e-10);

if aresame(str,'geo')
    [k2,l2,theta2,phi2]=elldiff(dt,kappa,lambda,theta,phi);
    v=ellrad(k2,l2,phi2);
elseif aresame(str,'ave')
    [k2,l2,theta2,phi2]=elldiff(dt,kappa,lambda,theta,phi);
    v=ellrad(k2,l2,phi2,'ave');
elseif aresame(str,'cir')
   omphi=frac(1,dt).*vdiff(unwrap(angle(rot(phi))),1,strend);
   omtheta=frac(1,dt).*vdiff(unwrap(angle(rot(theta))),1,strend);
   om=omphi+sign(lambda).*sqrt(1-lambda.^2).*omtheta;
   v=om;
   v=frac(om.*kappa,sqrt(1-lambda.^2));  %So that when I multiply by R it gives circulation, V_circ R = circ/(2pi)
elseif aresame(str,'kin')
    v=elldiff(dt,kappa,lambda,theta,phi);
elseif aresame(str,'azi')
    z=ellsig(kappa,lambda,theta,phi);
    Omega=frac(1,dt).*vdiff(unwrap(angle(z)),1,strend);
    %vsize(z,Omega,kappa,phi)
    v=Omega.*ellrad(kappa,lambda,phi,'instantaneous');
elseif aresame(str,'ins')
    [k2,l2,theta2,phi2]=elldiff(dt,kappa,lambda,theta,phi);
    v=ellrad(k2,l2,phi2,'instantaneous');
end

v=abs(v).*sign(lambda).*fact;    



