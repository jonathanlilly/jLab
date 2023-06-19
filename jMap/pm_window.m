function[varargout]=pm_window(varargin)
%PM_WINDOW  Remove data outside of time window for POLYMAP.
%
%   [XI,YI,TI]=PM_WINDOW(X,Y,T,MU), where the X, Y, and T are arrays of the
%   same size and MU is a temporal bandwith, returns indexed versions XI, 
%   YI, and TI corresponding to all times T such that ABS(T)<MU.
%  
%   [XI,YI,TI,Z1I,Z2I,...,ZNI]=PM_WINDOW(X,Y,T,Z1,Z2,...,ZN,MU) for any
%   number of additional input arguments also works.
%
%   Note that any empty input arrays are returned as empty.
%
%   PM_WINDOW is called internally by POLYMAP.  However, for large problems
%   it may be preferable to call it externally, as documented in POLYMAP.
%
%   See also POLYMAP.
%
%   Usage: [xi,yi,ti]=pm_window(x,y,t,mu);
%          [xi,yi,ti,zi]=pm_window(x,y,t,z,mu);
%          [xi,yi,ti,wi,zi]=pm_window(x,y,t,w,z,mu);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 
mu=varargin{end};
tdata=varargin{3};

%remove data points that are outide the temporal window
if ~isempty(tcell)
    index=abs(tdata./mu)<1;
    for i=1:length(nargin)-1
        if ~isempty(varargin{i})
            varargout{i}=varargin{i}(index);
        else
            varargout{i}=[];
        end
    end
end
 
