function[varargout]=arrayify(varargin)
%ARRAYIFY  Converts a set of scalars or arrays into column arrays.
%
%   [O1,O2,O3,...,ON]=ARRAYIFY(I1,I2,I3,...,IN) where each of the IN is
%   either a scalar or an array with M elements, converts each input 
%   variable into a length M column vector.  
%
%   The input variables that are M-element arrays are resized to column
%   vectors, while scalars are replicated to have M identical elements.
%
%   ARRAYIFY will return an error if any of the input variables that are 
%   not scalars have different numbers of elements from each other. 
%
%   ARRAYIFY(I1,I2,...IN); with no output arguments overwrites the original
%   input variables.
%
%   'arrayify --t' runs a test.
%
%   Usage: [o1,o2,o3]=arrayify(i1,i2,i3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2016 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    arrayify_test,return
end


M=zeros(length(varargin),1);
for i=1:length(varargin)
    varargin{i}=varargin{i}(:);
    M(i)=length(varargin{i});
end

if any(M>1)
    if any(M(M~=1)~=max(M))
        error('Input parameters must either be scalars or arrays having the same number of elements.')
    end
    M=max(M);
    for i=1:length(varargin)
        varargin{i}=varargin{i}+0*ones(M,1);
    end
    varargout=varargin;
    eval(to_overwrite(nargin));
else
    %Do nothing.  No need to overwrite in this case. 
    varargout=varargin;
end

function[]=arrayify_test

[o1,o2,o3]=arrayify(4,1:10,[1:10]');

reporttest('ARRAYIFY',aresame(o1,4+0*[1:10]')&aresame(o2,[1:10]')&aresame(o3,[1:10]'))

i1=4;
i2=1:10;
i3=[1:10]';

arrayify(i1,i2,i3);

reporttest('ARRAYIFY with overwriting',aresame(i1,4+0*[1:10]')&aresame(i2,[1:10]')&aresame(i3,[1:10]'))
