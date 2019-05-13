function[f]=fourier(varargin)
%FOURIER  Returns the Fourier frequencies for a given length time series.
%
%   F=FOURIER(M) returns the one-sided (or positive) Fourier frequencies 
%   for a time series of length M.   
%
%   F is a radian or angular frequency so that the Nyquist is at PI.
%
%   F=FOURIER(M,'two') instead returns the two-sided Fourier frequencies.
%   The default behavior is equivalent to F=FOURIER(M,'one').
%
%   F=FOURIER(DT,M) uses sample rate DT in calculating the frequencies, so
%   that the Nyquist will be at PI/DT.
% 
%   Note that the highest resolved frequency, MAX(F), differs for even or
%   odd M.  For even M, it is the Nyquist PI/DT, but for odd M the Nyquist
%   is not resolved and the highest resolved frequency is (M-1)/M * PI/DT.
%
%   For the one-sided option, F has length FLOOR(M/2)+1, or M/2+1 for even 
%   M, and (M+1)/2 for odd M.  For the two-sided option, F has length M.
%   __________________________________________________________________
%
%   Array input
%
%   F=FOURIER(M) also works if M an array instead of a scalar.  In this
%   case, F is a cell array with LENGTH(M) elements.
%
%   F=FOURIER(DT,M) with M being an array works provided DT is either a 
%   scalar or an array of the same length as M.   
%   __________________________________________________________________
%
%   See also CENTEREDTIMES.
%
%   Usage: f=fourier(M);
%          f=fourier(dt,M);
%          f=fourier(dt,M,'two');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2019 J.M. Lilly --- type 'help jlab_license' for details
 
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
