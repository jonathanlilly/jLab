function[varargout]=cellget(varargin)
%CELLGET  Indexes a cell array of numerical arrays by ID number.
%
%   [Y1,Y2,...,YN]=CELLGET(X1,X2,...,XN,ID,IDO) where the N input arrays
%   are cell arrays of column vectors, returns those elements of the XN for
%   which the identification number ID corresponds to requested list IDO.
%
%   ID is a numerical array of the same length as the input cell arrays.
%   IDO may either be a scalar or an array. 
%
%   If IDO is an array, then the output fields Y1,Y2,...,YN are cell arrays
%   with the same length as IDO, containing the corresponding elements of 
%   the input fields X1,X2,...,XN.  If a particular entry of IDO is not 
%   found within ID, that element of Y1,Y2,...,YN is left empty.  
%
%   If IDO is an array, then Y1,Y2,...,YN are numerical arrays containing 
%   only the requested element of the input field.  The YN will be empty if 
%   the requested identification number is not found within ID.
%
%   CELLGET(X1,X2,...,XN,ID,IDO); with no output arguments overwrites the 
%   named input variables X1,X2,...,XN, and also overwrites ID with IDO.
%   _______________________________________________________________________
%
%   Use with Lagrangian data
%
%   CELLGET is primarily used to extract float or drifter trajectories by 
%   identification number, in particular, with FLOATS.MAT and DRIFTERS.MAT.  
%
%   For example, with DRIFTER.MAT in memory, USE DRIFTERS followed by 
%   CELLGET(NUM,LAT,LON,ID,44000) overwrites cell arrays NUM, LAT, and 
%   LON as numerical arrays containing the trajectory of drifter ID#44000. 
%   This particular drifter has come to be known as 'drifter Betty'.
%   _______________________________________________________________________
%
%   Usage: [y1,y2,y3]=cellget(x1,x2,x3,id,ido);
%          cellget(x1,x2,x3,id,ido);
%          cellget(num,lat,lon,id,ido);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
ido=varargin{end};
id=varargin{end-1};
varargin=varargin(1:end-2);


if length(ido)==1
    ii=find(id==ido);
    for i=1:length(varargin)
        if ~isempty(ii)
            varargout{i}=varargin{i}{ii};
        else
            varargout{i}=[];
        end
    end
elseif length(ido)>1
    for j=1:length(ido)
        ii=find(id==ido(j));
        for i=1:length(varargin)
            if ~isempty(ii)
                varargout{i}{j,1}=varargin{i}{ii};
            else
                varargout{i}{j,1}=[];
            end
        end
    end
end

varargout{end+1}=ido;
eval(to_overwrite(length(varargout)));