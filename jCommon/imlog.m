function[y]=imlog(x)
%IMLOG  Imaginary part of log: IMLOG(X)=UNWRAP(IMAG(LOG(X)))
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details  

y=unwrap(imag(log(x)));
