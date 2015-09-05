function[varargout]=trajwrap(varargin)
%TRAJWRAP  Wraps Lagrangian trajectories to fit within a periodic domain.
%
%   CX=TRAJWRAP(CX,L) wraps the complex-valued Lagrangian trajectory CX to 
%   fit within the doubly periodic domain with edges at +/- L.  
%
%   CX specifies the particle position in complex-valued form CX=X+iY.
%
%   CX may be of any dimension, but time is assumed to increase along rows.
%   The output array will be the same size as the input array.  
%
%   Set L=[LX LY] to specify different periodicities in the X and Y 
%   directions. Use L=[LX NaN] for only zonal periodicity, or L=[NaN Ly] 
%   for only meridional periodicity.
%
%   CX=TRAJWRAP(CX,L,'nans') optionally sets the discontinuties resulting
%   from wrapping to values of NaN, for plotting purposes.  Note that the 
%   output is in this case still the same size as the input. 
%
%   [CX1,CX2,...,CXN]=TRAJWRAP(CX1,CX2,...CXN,L) also works for multiple
%   arrays of trajectories with the same periodicity.
%
%   TRAJWRAP(CX1,CX2,...CXN,L); with no output arguments overwrites the 
%   original named input variables.
%
%   See also TRAJUNWRAP. 
%
%   'trajwrap --t' runs a test.
%   'trajwrap --f' generates a sample figure.
%
%   Usage: cx=trajwrap(cx,L);
%          cx=trajwrap(cx,[Lx Ly],'nan');
%          cx=trajwrap(cx,[Lx Ly]);
%          [cx1,cx2]=trajwrap(cx1,cx2,[Lx Ly]);
%          trajwrap(cx1,cx2,[Lx Ly]);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}(1:3),'--t')
    trajwrap_test;return
elseif strcmpi(varargin{1}(1:3),'--f')
    type makefigs_trajwrap
    makefigs_trajwrap;
    return
end

str='nobreaks';
na=nargin;
if ischar(varargin{end})
    str=varargin{end};
    na=na-1;
    varargin=varargin(1:end-1);
end

L=varargin{end};
if length(L)==1
    Lx=L;
    Ly=L;
else
   Lx=L(1);
   Ly=L(2);
end

for i=1:na-1
    z=varargin{i};
    if iscell(z)
        for j=1:length(z)
            varargout{i}{j}=trajwrap_one(Lx,Ly,z{j},str);
        end
    else
       varargout{i}=trajwrap_one(Lx,Ly,z,str);
    end
end

eval(to_overwrite(na-1));         

function[zw]=trajwrap_one(Lx,Ly,z,str)
x=real(z);
y=imag(z);
%Lx,Ly
if ~isnan(Lx)
    x=frac(Lx,pi)*angle(rot(frac(pi,Lx)*real(z)));
end
if ~isnan(Ly)
    y=frac(Ly,pi)*angle(rot(frac(pi,Ly)*imag(z)));
end
zw=x+sqrt(-1)*y;
if ~isempty(strfind(str,'nan'))
    zwshift=vshift(zw,1,1);
    bool=(abs(real(zw-zwshift))>Lx/2)|(abs(imag(zw-zwshift))>Ly/2);
    bool(1,:,:,:,:,:,:)=true;
    bool(end,:,:,:,:,:,:)=true;
    zw(bool)=nan+sqrt(-1)*nan;
end

function[]=trajwrap_test

load vortex
use vortex.drifters
cx=x+1i*y;
cx2=trajwrap(trajunwrap(cx,[L*pi,nan]),[L*pi,nan]);
reporttest('TRAJWRAP inverts TRAJUNWRAP for VORTEX dataset',aresame(cx2,cx,1e-10))




