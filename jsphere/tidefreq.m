function[f,name]=tidefreq(str)
%TIDEFREQ  Frequencies of the eight major tidal components.
%
%   [F,NAME]=TIDEFREQ returns the frequencies of the eight major tidal 
%   components together with a string matrix containing their names.  In
%   order of increasing frequency, these are   
%    
%        Mf O1 P1 K1 N2 M2 S2 K2.
%
%   The units are given in *radians* per hour. See Gill (1982) page 335.
%  
%   F=TIDEFREQ(STR) returns the frequency for the component named STR, 
%   if STR is a string.  F=TIDEFREQ(N) where N is a number also works.
%  
%   Usage:   [f,name]=tidefreq;
%            f=tidefreq(str);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2015 J.M. Lilly --- type 'help jlab_license' for details    

name=('MFO1P1K1N2M2S2K2');
  
%Gill only gives first two decimal points
%I got the more exact numbers from
%http://oceanworld.tamu.edu/resources/ocng_textbook/chapter17/chapter17_04.htm

p(1)=327.85;
p(2)=25.8194;
p(3)=24.0659;
p(4)=23.9344;
p(5)=12.6584;
p(6)=12.4206;
p(7)=12.0000;
p(8)=11.9673;

f=1./p;
if nargin==1
  if ~ischar(str)
    f=f(str);
  else
    str=upper(str);
    ii=strfind(name,str);
    f=f((ii-1)/2+1); 
  end
  clear name
else
  name=reshape(name,2,8)';
end

f=f*2*pi;