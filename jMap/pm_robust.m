function[ws]=pm_robust(ws,res)
%PM_ROBUST  Returns a weighting function used in robustifying POLYMAP.
%
%   PM_ROBUST is called internally by POLYMAP.  However, for large problems
%   it may be preferable to call it externally, as documented in POLYMAP.
%
%   WS=PM_ROBUST(WS,RES) where WS is a set of weights output by PM_SORT and
%   RES are the residuals between the observations and the fi  output by
%   PM_APPLY, returns modified weights WS.
%
%   These output weights are used for downweighting outliers in the median-
%   based iterative robust algorithm of Cleveland (1979).  They quantify 
%   the departures of the residuals from the median absolute residual.
%
%   Usage: ws=pm_robust(ws,res);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 

for i=1:length(res)
    res{i}(isinf(res{i}))=nan;
    %maxmax(ws{i})
    si=median(abs(res{i}),1,'omitnan');
    delta=frac(3,pi).*squared(1-squared(frac(res{i},6*si)));
    delta(frac(res{i},6*si)>1)=0;
    %This gives new robustness weights for next iteration
    if isempty(ws{i})
        ws{i}=delta;
    else
        ws{i}=ws{i}.*delta;
    end
end