function z = jdawson(x,n)
%JDAWSON The Dawson function and its derivatives. [With P.J. Acklam]
%
%   Y = JDAWSON(X) is Dawson's integral for each element of X, with X real.
%   
%   Dawson's integral F(x) is defined as
%
%     F(x) = exp(-x^2) * integral from 0 to x of exp(t^2) dt
%
%   The Dawson function is also the Hilbert transform of the Gaussian,
%   see for example Lilly and Olhede (2009).
%
%   This function was written almost entirely by P.J. Acklam and is 
%   redistributed with permission.
%   __________________________________________________________________
%
%   Derivatives
%
%   YN=JDAWSON(X,N) optionally returns the Nth derivatives of the Dawson 
%   function at positions X. 
%
%   A form for the derivatives of the JDAWSON function is given by Lilly 
%   and Olhede (2008b), and involves Hermite polynomials.
%   __________________________________________________________________
%  
%   See also ERF, ERFC, HERMPOLY.  
%
%   'jdawson --t' runs a test.
%   'jdawson --f' generates a sample figure.
%
%   Usage: y=jdawson(x);
%          yn=jdawson(x,n);
%   __________________________________________________________________
%   (C) 2004--2016  P.J. Acklam and J.M. Lilly 
%                             --- type 'help jlab_license' for details

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-02-09 16:21:09 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

%   This is a MATLAB implementation of a FORTRAN routine which is available
%   at <URL: http://www.netlib.org/fn/ddaws.f >.

if strcmpi(x, '--t')
    jdawson_test,return
end
if strcmpi(x, '--f')
    type makefigs_jdawson
    makefigs_jdawson;
    return
end

if nargin==1
    z=jdawsonfun(x);
else
    z=jdawsonderiv(x,n);
end

function[y]=jdawsonderiv(x,n)    

herm=hermpoly(x(:),n);
iherm=hermpoly(sqrt(-1)*x(:),n);

for k=1:n+1
    hermcell{k}=reshape(herm(:,k),size(x));
    ihermcell{k}=reshape(iherm(:,k),size(x));
    %real(sqrt(-1).^(k-1).*reshape(iherm(:,k),size(x)));
end
y=zeros(size(x));
for k=1:n
     y=y+choose(n,k).*hermcell{n-k+1}.*((sqrt(-1).^(k-1)).*ihermcell{k});
     %                      H_(n-k}                               H{k-1}
end
y=((-1).^n).*(hermcell{n+1}.*jdawson(x)-y);


function[z]=jdawsonfun(x)    
% Check array class.
if ~isa(x, 'double')
    error('Input must be a double array.');
end

% Check array values.
if ~isreal(x)
    error('Input must be real.');
end

%
% may 1980 edition.  w. fullerton c3, los alamos scientific lab.
%

% series for daw        on the interval  0.          to  1.00000e+00
%                                        with weighted error   8.95e-32
%                                         log weighted error  31.05
%                               significant figures required  30.41
%                                    decimal places required  31.71
%
dawcs = zeros(21, 1);
dawcs( 1) = -.6351734375145949201065127736293d-2;
dawcs( 2) = -.2294071479677386939899824125866d+0;
dawcs( 3) = +.2213050093908476441683979161786d-1;
dawcs( 4) = -.1549265453892985046743057753375d-2;
dawcs( 5) = +.8497327715684917456777542948066d-4;
dawcs( 6) = -.3828266270972014924994099521309d-5;
dawcs( 7) = +.1462854806250163197757148949539d-6;
dawcs( 8) = -.4851982381825991798846715425114d-8;
dawcs( 9) = +.1421463577759139790347568183304d-9;
dawcs(10) = -.3728836087920596525335493054088d-11;
dawcs(11) = +.8854942961778203370194565231369d-13;
dawcs(12) = -.1920757131350206355421648417493d-14;
dawcs(13) = +.3834325867246327588241074439253d-16;
dawcs(14) = -.7089154168175881633584099327999d-18;
dawcs(15) = +.1220552135889457674416901120000d-19;
dawcs(16) = -.1966204826605348760299451733333d-21;
dawcs(17) = +.2975845541376597189113173333333d-23;
dawcs(18) = -.4247069514800596951039999999999d-25;
dawcs(19) = +.5734270767391742798506666666666d-27;
dawcs(20) = -.7345836823178450261333333333333d-29;
dawcs(21) = +.8951937667516552533333333333333d-31;
%
% series for daw2       on the interval  0.          to  1.60000e+01
%                                        with weighted error   1.61e-32
%                                         log weighted error  31.79
%                               significant figures required  31.40
%                                    decimal places required  32.62
%
daw2cs = zeros(45, 1);
daw2cs( 1) = -.56886544105215527114160533733674d-1;
daw2cs( 2) = -.31811346996168131279322878048822d+0;
daw2cs( 3) = +.20873845413642236789741580198858d+0;
daw2cs( 4) = -.12475409913779131214073498314784d+0;
daw2cs( 5) = +.67869305186676777092847516423676d-1;
daw2cs( 6) = -.33659144895270939503068230966587d-1;
daw2cs( 7) = +.15260781271987971743682460381640d-1;
daw2cs( 8) = -.63483709625962148230586094788535d-2;
daw2cs( 9) = +.24326740920748520596865966109343d-2;
daw2cs(10) = -.86219541491065032038526983549637d-3;
daw2cs(11) = +.28376573336321625302857636538295d-3;
daw2cs(12) = -.87057549874170423699396581464335d-4;
daw2cs(13) = +.24986849985481658331800044137276d-4;
daw2cs(14) = -.67319286764160294344603050339520d-5;
daw2cs(15) = +.17078578785573543710504524047844d-5;
daw2cs(16) = -.40917551226475381271896592490038d-6;
daw2cs(17) = +.92828292216755773260751785312273d-7;
daw2cs(18) = -.19991403610147617829845096332198d-7;
daw2cs(19) = +.40963490644082195241210487868917d-8;
daw2cs(20) = -.80032409540993168075706781753561d-9;
daw2cs(21) = +.14938503128761465059143225550110d-9;
daw2cs(22) = -.26687999885622329284924651063339d-10;
daw2cs(23) = +.45712216985159458151405617724103d-11;
daw2cs(24) = -.75187305222043565872243727326771d-12;
daw2cs(25) = +.11893100052629681879029828987302d-12;
daw2cs(26) = -.18116907933852346973490318263084d-13;
daw2cs(27) = +.26611733684358969193001612199626d-14;
daw2cs(28) = -.37738863052129419795444109905930d-15;
daw2cs(29) = +.51727953789087172679680082229329d-16;
daw2cs(30) = -.68603684084077500979419564670102d-17;
daw2cs(31) = +.88123751354161071806469337321745d-18;
daw2cs(32) = -.10974248249996606292106299624652d-18;
daw2cs(33) = +.13261199326367178513595545891635d-19;
daw2cs(34) = -.15562732768137380785488776571562d-20;
daw2cs(35) = +.17751425583655720607833415570773d-21;
daw2cs(36) = -.19695006967006578384953608765439d-22;
daw2cs(37) = +.21270074896998699661924010120533d-23;
daw2cs(38) = -.22375398124627973794182113962666d-24;
daw2cs(39) = +.22942768578582348946971383125333d-25;
daw2cs(40) = -.22943788846552928693329592319999d-26;
daw2cs(41) = +.22391702100592453618342297600000d-27;
daw2cs(42) = -.21338230616608897703678225066666d-28;
daw2cs(43) = +.19866196585123531518028458666666d-29;
daw2cs(44) = -.18079295866694391771955199999999d-30;
daw2cs(45) = +.16090686015283030305450666666666d-31;
%
% series for dawa       on the interval  0.          to  6.25000e-02
%                                        with weighted error   1.97e-32
%                                         log weighted error  31.71
%                               significant figures required  29.79
%                                    decimal places required  32.64
%
dawacs = zeros(75, 1);
dawacs( 1) = +.1690485637765703755422637438849d-1;
dawacs( 2) = +.8683252278406957990536107850768d-2;
dawacs( 3) = +.2424864042417715453277703459889d-3;
dawacs( 4) = +.1261182399572690001651949240377d-4;
dawacs( 5) = +.1066453314636176955705691125906d-5;
dawacs( 6) = +.1358159794790727611348424505728d-6;
dawacs( 7) = +.2171042356577298398904312744743d-7;
dawacs( 8) = +.2867010501805295270343676804813d-8;
dawacs( 9) = -.1901336393035820112282492378024d-9;
dawacs(10) = -.3097780484395201125532065774268d-9;
dawacs(11) = -.1029414876057509247398132286413d-9;
dawacs(12) = -.6260356459459576150417587283121d-11;
dawacs(13) = +.8563132497446451216262303166276d-11;
dawacs(14) = +.3033045148075659292976266276257d-11;
dawacs(15) = -.2523618306809291372630886938826d-12;
dawacs(16) = -.4210604795440664513175461934510d-12;
dawacs(17) = -.4431140826646238312143429452036d-13;
dawacs(18) = +.4911210272841205205940037065117d-13;
dawacs(19) = +.1235856242283903407076477954739d-13;
dawacs(20) = -.5788733199016569246955765071069d-14;
dawacs(21) = -.2282723294807358620978183957030d-14;
dawacs(22) = +.7637149411014126476312362917590d-15;
dawacs(23) = +.3851546883566811728777594002095d-15;
dawacs(24) = -.1199932056928290592803237283045d-15;
dawacs(25) = -.6313439150094572347334270285250d-16;
dawacs(26) = +.2239559965972975375254912790237d-16;
dawacs(27) = +.9987925830076495995132891200749d-17;
dawacs(28) = -.4681068274322495334536246507252d-17;
dawacs(29) = -.1436303644349721337241628751534d-17;
dawacs(30) = +.1020822731410541112977908032130d-17;
dawacs(31) = +.1538908873136092072837389822372d-18;
dawacs(32) = -.2189157877645793888894790926056d-18;
dawacs(33) = +.2156879197938651750392359152517d-20;
dawacs(34) = +.4370219827442449851134792557395d-19;
dawacs(35) = -.8234581460977207241098927905177d-20;
dawacs(36) = -.7498648721256466222903202835420d-20;
dawacs(37) = +.3282536720735671610957612930039d-20;
dawacs(38) = +.8858064309503921116076561515151d-21;
dawacs(39) = -.9185087111727002988094460531485d-21;
dawacs(40) = +.2978962223788748988314166045791d-22;
dawacs(41) = +.1972132136618471883159505468041d-21;
dawacs(42) = -.5974775596362906638089584995117d-22;
dawacs(43) = -.2834410031503850965443825182441d-22;
dawacs(44) = +.2209560791131554514777150489012d-22;
dawacs(45) = -.5439955741897144300079480307711d-25;
dawacs(46) = -.5213549243294848668017136696470d-23;
dawacs(47) = +.1702350556813114199065671499076d-23;
dawacs(48) = +.6917400860836148343022185660197d-24;
dawacs(49) = -.6540941793002752512239445125802d-24;
dawacs(50) = +.6093576580439328960371824654636d-25;
dawacs(51) = +.1408070432905187461501945080272d-24;
dawacs(52) = -.6785886121054846331167674943755d-25;
dawacs(53) = -.9799732036214295711741583102225d-26;
dawacs(54) = +.2121244903099041332598960939160d-25;
dawacs(55) = -.5954455022548790938238802154487d-26;
dawacs(56) = -.3093088861875470177838847232049d-26;
dawacs(57) = +.2854389216344524682400691986104d-26;
dawacs(58) = -.3951289447379305566023477271811d-27;
dawacs(59) = -.5906000648607628478116840894453d-27;
dawacs(60) = +.3670236964668687003647889980609d-27;
dawacs(61) = -.4839958238042276256598303038941d-29;
dawacs(62) = -.9799265984210443869597404017022d-28;
dawacs(63) = +.4684773732612130606158908804300d-28;
dawacs(64) = +.5030877696993461051647667603155d-29;
dawacs(65) = -.1547395051706028239247552068295d-28;
dawacs(66) = +.6112180185086419243976005662714d-29;
dawacs(67) = +.1357913399124811650343602736158d-29;
dawacs(68) = -.2417687752768673088385304299044d-29;
dawacs(69) = +.8369074582074298945292887587291d-30;
dawacs(70) = +.2665413042788979165838319401566d-30;
dawacs(71) = -.3811653692354890336935691003712d-30;
dawacs(72) = +.1230054721884951464371706872585d-30;
dawacs(73) = +.4622506399041493508805536929983d-31;
dawacs(74) = -.6120087296881677722911435593001d-31;
dawacs(75) = +.1966024640193164686956230217896d-31;
%
eps = d1mach(3);
ntdaw  = initds(dawcs,  21, 0.1*eps);
ntdaw2 = initds(daw2cs, 45, 0.1*eps);
ntdawa = initds(dawacs, 75, 0.1*eps);
%
xsml = sqrt(1.5*eps);
xbig = sqrt(0.5/eps);
xmax = exp(min(-log(2.d0*d1mach(1)), log(d1mach(2))) - 0.01d0);
%
z = zeros(size(x));
%
y = abs(x);
%
k = (y <= xsml);
if any(k(:))
    z(k) = x(k);
end
%
k = (xsml < y) & (y <= 1.0d0);
if any(k(:))
    z(k) = x(k) .* (.75d0 + dcsevl(2.d0*y(k).^2-1.d0, dawcs, ntdaw));
end
%
k = (1.0d0 < y) & (y <= 4.d0);
if any(k(:))
    z(k) = x(k) .* (.25d0 + dcsevl(.125d0*y(k).^2-1.d0, daw2cs, ntdaw2));
end
%
k = (4.d0 < y) & (y <= xbig);
if any(k(:))
    z(k) = (0.5d0 + dcsevl(32.d0./y(k).^2-1.d0, dawacs, ntdawa)) ./ x(k);
end
%
k = (xbig < y) & (y <= xmax);
if any(k(:))
    z(k) = 0.5d0./x(k);
end
%
k = (xmax < y);
if any(k(:))
    z(k) = 0.0d0;
end
%
k = isnan(x);
if any(k(:))
    z(k) = NaN;
end

function d = d1mach(i)

%   This is a MATLAB implementation of a FORTRAN routine which is available
%   at <URL: http://www.netlib.org/slatec/src/d1mach.f >.

%   D1MACH can be used to obtain machine-dependent parameters for the local
%   machine environment.  It is a function subprogram with one (input)
%   argument, and can be referenced as follows:
%
%        D = D1MACH(I)
%
%   where I=1,...,5.  The (output) value of D above is determined by the
%   (input) value of I.  The results for various values of I are discussed
%   below.
%
%        D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
%        D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
%        D1MACH( 3) = B**(-T), the smallest relative spacing.
%        D1MACH( 4) = B**(1-T), the largest relative spacing.
%        D1MACH( 5) = LOG10(B)
%
%   Assume double precision numbers are represented in the T-digit, base-B
%   form
%
%              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
%
%   where 0 <= X(I) < B for I=1,...,T, 0 < X(1), and EMIN <= E <= EMAX.
%
%   The values of B, T, EMIN and EMAX are provided as follows:
%
%        B    = the base.
%        T    = the number of base-B digits.
%        EMIN = the smallest exponent E.
%        EMAX = the largest exponent E.

dmach = zeros(5, 1);
dmach(1) = realmin;
dmach(2) = realmax;
dmach(3) = eps / 2;
dmach(4) = eps;
dmach(5) = .301029995663981195213738894724493027e+000;   % log(2)/log(10)

d = dmach(i);

function initds = initds(dos, nos, eta)

%   This is a MATLAB implementation of a FORTRAN routine which is available
%   at <URL: http://www.netlib.org/fn/initds.f >.

% june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
%
% initialize the double precision orthogonal series dos so that initds
% is the number of terms needed to insure the error is no larger than
% eta.  ordinarily eta will be chosen to be one-tenth machine precision.
%
%             input arguments --
% dos    dble prec array of nos coefficients in an orthogonal series.
% nos    number of coefficients in dos.
% eta    requested accuracy of series.

if (nos < 1)
    error('number of coefficients < 1');
end

err = 0.;
for ii = 1:nos
    i = nos + 1 - ii;
    err = err + abs(sngl(dos(i)));
    if (err > eta)
        break;
    end
end

if (i == nos)
    error('eta may be too small');
end
initds = i;

function y = sngl(x)

%   This is a replacement function for the FORTRAN intrinsic function
%   `sngl'.  MATLAB does not do arithmetic on single precision numbers, so
%   convert to single and then back again to double.

y = double(single(x));

function dcsevl = dcsevl(x, a, n)

%   This is a MATLAB implementation of a FORTRAN routine which is available
%   at <URL: http://www.netlib.org/fn/dcsevl.f >.

% evaluate the n-term chebyshev series a at x.  adapted from
% r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
%
%             input arguments --
% x      dble prec value at which the series is to be evaluated.
% a      dble prec array of n terms of a chebyshev series.  in eval-
%        uating a, only half the first coef is summed.
% n      number of terms in array a.

if (n < 1)
    error('number of terms <= 0');
end
if (n > 1000)
    error('number of terms > 1000');
end
if (any(x < (-1.1d0)) || any(x > 1.1d0))
    error('x outside (-1,+1)');
end

twox = 2.0d0*x;
b1 = 0.d0;
b0 = 0.d0;
for i = 1:n
    b2 = b1;
    b1 = b0;
    ni = n - i + 1;
    b0 = twox.*b1 - b2 + a(ni);
end

dcsevl = 0.5d0 * (b0-b2);
        
function[]=jdawson_test

dt=0.1;
t=(-150:dt:150)';
g=exp(-t.^2);
psi1=frac(2,sqrt(pi))*jdawson(t);
%plot(imag(anatrans(g))-imag(hilbert(g))),hold on

%aresame(hilbert(g),imag(anatrans(g)),1e-10)
psi2=imag(anatrans(g));

err=vsum(abs(psi1-psi2).^2,1)./vsum(abs(psi1).^2,1);
%plot(t,[psi1 psi2])
reporttest('JDAWSON equals Hilbert transform of Gaussian',err<1e-3)

dt=0.01;
t=(-15:dt:15)';
dk1=jdawson(t);
dk1([1 end],:)=nan;  %Avoid differentiation errors at edges
for k=1:5
    dk2=jdawson(t,k);
    dk1=vdiff(dk1,1)./dt;
    err=vsum(abs(dk1-dk2).^2,1)./vsum(abs(dk1).^2,1);
    reporttest(['JDAWSON derivatives match numerical differentiation for n=' int2str(k)],err<1e-3)
end


  