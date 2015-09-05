function[varargout]=cum2mom(varargin)
%CUM2MOM  Convert cumulants to moments.
%
%   [M0,M1,...MN]=CUM2MOM(K0,K1,...KN) converts the first N cumulants 
%   K0,K1,...KN into the first N moments M0,M1,...MN.
%
%   The KN and MN are all scalars or arrays of the same size.
%
%   Note for a probability density function, M0=1 and K0=0.
%
%   MCELL=CUM2MOM(KCELL), where KCELL is a cell array whose (N+1)th
%   is the Nth cumulant, returns a similar cell array of moments.
%
%   CUM2MOM is inverted by MOM2CUM.  
%
%   See also MOM2CUM, BELLPOLY.
%  
%   Usage: [m0,m1,m2]=cum2mom(k0,k1,k2);
%          mcell=cum2mom(kcell);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details

%Tests for CUM2MOM are run by MOM2CUM.
if strcmpi(varargin{1}, '--t'),return,end

if nargin==1&&iscell(varargin{1})
    cum=varargin{1};
else
    cum=varargin;
end
mom=cell(size(cum));
mom{1}=exp(cum{1});

cum=cum(2:end);
mom(2:end)=bellpoly(cum);

for i=2:length(mom)
    mom{i}=mom{i}.*mom{1};
end

if nargin==1&&iscell(varargin{1})
    varargout{1}=mom;
else
    varargout=mom;
end

