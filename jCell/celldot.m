function[p]=celldot(x,y,str)
%CELLDOT  Dot product for arrays of column vectors.
%
%   P=CELLDOT(X,Y), with X and Y being cell arrays of column vectors,
%   returns a column vector P of length LENGTH(X), such that the nth 
%   element of P is given by P(n)=SUM(X{n}.*CONJ(Y{n})).
%
%   Note that this dot product is defined to conjugate the second input 
%   argument, unlike Matlab's DOT which conjugates the first argument.
%
%   P=CELLDOT(X,Y,'norm') normalizes the magnitude of dot product to be
%   between zero and one.  This is done by dividing the nth element by
%   SQRT(SUM(ABS(X{n}).^2).*SUM(ABS(Y{n}).^2)).
%
%   CELLDOT uses 'omitnan' to exclude any NaNs from the summation.
%
%   Usage: p=celldot(x,y);
%          p=celldot(x,y,'norm');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details
 
% if strcmp(x, '--t')
%     celldot_test,return
% end

if nargin ==2 
    str='nonorm';
end

%   'celldot --t' runs a test.
%if ~iscell(x) 
%    p=celldot_one(x,y);
%else
    p=zeros(length(x),1);
    for i=1:length(x)
        p(i)=celldot_one(x{i},y{i});
    end
%end

if strcmp(str(1:3),'nor')
    %if ~iscell(x)
    %    p=p./sqrt(celldot_one(x,x).*celldot_one(y,y));
    %else
        for i=1:length(x)
            p(i)=p(i)./sqrt(celldot_one(x{i},x{i}).*celldot_one(y{i},y{i}));
        end
    %end
end


function[p]=celldot_one(x,y)

p=sum(x.*conj(y),'omitnan');

%function[]=celldot_test
 
%reporttest('CELLDOT',aresame())
