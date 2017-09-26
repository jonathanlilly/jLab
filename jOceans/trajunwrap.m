function[varargout]=trajunwrap(varargin)
%TRAJUNWRAP  Unwraps Lagrangian trajectories from a periodic domain.
%
%   CX=TRAJUNWRAP(CX,L) unwraps the complex-valued Lagrangian trajectory CX 
%   from within a doubly periodic domain having edges at +/- L.  
%
%   CX specifies the partile position in complex-valued form CX=X+iY.
%
%   CX may be of any dimension, but time is assumed to increase along rows.
%   The output array will be the same size as the input array.  
%
%   Set L=[LX LY] to specify different periodicities in the X and Y 
%   directions. Use L=[LX NaN] for only zonal periodicity, or L=[NaN Ly] 
%   for only meridional periodicity.
%
%   [CX1,CX2,...,CXN]=TRAJUNWRAP(CX1,CX2,...CXN,L) also works for multiple
%   arrays of trajectories with the same periodicity.
%
%   TRAJUNWRAP(CX1,CX2,...CXN,L); with no output arguments overwrites the 
%   original named input variables.
%
%   See also TRAJWRAP.
%
%   Usage: cx=trajunwrap(cx,L);
%          cx=trajunwrap(cx,[Lx Ly],'nan');
%          cx=trajunwrap(cx,[Lx Ly]);
%          [cx1,cx2]=trajunwrap(cx1,cx2,[Lx Ly]);
%          trajunwrap(cx1,cx2,[Lx Ly]);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 
na=nargin;
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
            varargout{i}{j}=trajunwrap_one(Lx,Ly,z{j});
        end
    else
       varargout{i}=trajunwrap_one(Lx,Ly,z);
    end
end

eval(to_overwrite(na-1));         

function[zw]=trajunwrap_one(Lx,Ly,z)
x=real(z);
y=imag(z);
if ~isnan(Lx)
    x=frac(Lx,pi)*unwrap(frac(pi,Lx)*real(z),[],1);
end
if ~isnan(Ly)
    y=frac(Ly,pi)*unwrap(frac(pi,Lx)*imag(z),[],1);
end

zw=x+sqrt(-1)*y;

% if ~isempty(strfind(str,'nan'))
%     zwshift=vshift(zw,1,1);
%     bool=(abs(real(zw-zwshift))>Lx/2)|(abs(imag(zw-zwshift))>Ly/2);
%     bool(1,:,:,:,:,:,:)=true;
%     bool(end,:,:,:,:,:,:)=true;
%     zw(bool)=nan;
% end


