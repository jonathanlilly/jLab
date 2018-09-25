function[f]=fourier(varargin)
%FOURIER  Returns the Fourier frequencies for a given length time series.
%
%   F=FOURIER(N) returns the one-sided (or positive) Fourier frequencies 
%   for a time series of length N.   
%
%   F is a radian or angular frequency so that the Nyquist is at PI.
%
%   F=FOURIER(N,'two') instead returns the two-sided Fourier frequencies.
%   The default behavior is equivalent to F=FOURIER(N,'one').
%
%   F=FOURIER(DT,N) uses sample rate DT in calculating the frequencies, so
%   that the Nyquist will be at PI/DT.
% 
%   Note that the highest resolved frequency, MAX(F), differs for even or
%   odd N.  For even N, it is the Nyquist PI/DT, but for odd N the Nyquist
%   is not resolved and the highest resolved frequency is (N-1)/N * PI/DT.
%
%   For the one-sided option, F has length FLOOR(N/2)+1, or N/2+1 for even 
%   N, and (N+1)/2 for odd N.  For the two-sided option, F has length N.
%   __________________________________________________________________
%
%   Array input
%
%   F=FOURIER(N) also works if N an array instead of a scalar.  In this
%   case, F is a cell array with LENGTH(N) elements.
%
%   F=FOURIER(DT,N) with N being an array works provided DT is either a 
%   scalar or an array of the same length as N.   
%   __________________________________________________________________
%
%   See also CENTEREDTIMES.
%
%   Usage: f=fourier(N);
%          f=fourier(dt,N);
%          f=fourier(dt,N,'two');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2018 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    fourier_test,return
end

str='one';
if ischar(varargin{end})
    str=varargin{end}(1:3);
    varargin=varargin(1:end-1);
end

if length(varargin)==1
    N=varargin{1};
    dt=1;
else
    dt=varargin{1};
    N=varargin{2};
end

if length(N)==1
    f=fourier_one(dt,N,str);
else
    if length(dt)==1
        dt=dt+zeros(size(N));
    elseif length(dt)~=length(N)
        error('If N is an array, DT must either be a scalar or an array of the same length.')
    end
    for i=1:length(N)
        f{i,1}=fourier_one(dt,N(i),str);
    end
end

function[f]=fourier_one(dt,N,str)

if strcmpi(str,'one')
    f=2*pi*(0:floor(N/2))'./N;
else
    f=2*pi*(0:N-1)'./N;
    index=find(f>pi);
    f(index)=f(index)-2*pi;
end

% if iseven(N)
%     f=(0:1./N:1/2)';
% elseif isodd(N)
%     f=(0:1./N:1/2*frac(N-1,N))';
% end

f=f./dt;

function[]=fourier_test

reporttest('FOURIER is length N/2+1 for even N',length(fourier(10))==6);
reporttest('FOURIER is length (N+1)/2 for odd N',length(fourier(11))==6);

%Even and odd length frequencies behave differently, which you can see with 
%abs(fft(cos(pi*[1:10]'))).^2
%abs(fft(cos(pi*[1:11]'))).^2
%abs(fft(cos((10/11)*pi*[1:11]'))).^2
