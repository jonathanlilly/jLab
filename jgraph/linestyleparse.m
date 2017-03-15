function[widthcell,stylecell,colorcell]=linestyleparse(str)
%LINESTYLEPARSE  Parses the input string to LINESTYLE.
%
%   [WIDTH,STYLE,COLOR]=LINESTYLEPARSE(STR) converts the linestyle string 
%   STR into specifications for line style, width, and color.  It is a 
%   low-level function called by LINESTYLE and also TOPOPLOT.
%
%   COLORS=LINESTYLEPARSE; with no input arguments just outputs the 
%   definitions of colors used by LINESTYLE.
%
%   Usage: [width,style,color]=linestyleparse(str);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2016 J.M. Lilly --- type 'help jlab_license' for details
 
%Parse the linestyle input string for cell specifications

defaultwidth=1;
defaultcolor='k';
defaultstyle='-';
colors=linestyle_colors;

if nargin==0
    widthcell=colors;return
end

bool=true(size(str));

%/********************************************************
%linewidth
index=find(real(str)>47&real(str)<58); %find str is 0-9
if isempty(index)
   widthcell=defaultwidth;
else
   width=str2num(str(index(1):index(end)));  %encompasss '.'
   if isempty(width)
     error('Problem with linestyle string.')
   else
     widthcell=width;
     
     %remove these entries from string
     bool(index(1):index(end))=0;     
     str=str(bool);
     bool=true(size(str));
   end
end
%\********************************************************


%/********************************************************
%linecolor
index=find(real(str)>64&real(str)<123&~...
      (real(str)==real('s')|...
       real(str)==real('d')|...
       real(str)==real('v')|...
       real(str)==real('p')|...
       real(str)==real('h')|...
       real(str)==real('o')|...
       real(str)==real('x')));     %find str is a-Z not  s d v p h o x 
if isempty(index)
   colorcell=defaultcolor;
else
   color=str(index); 
   if length(color)>1
     error('More than one color specified.')
   else
     colorcell=getfield(colors,color);
     
     %remove these entries from string
     bool(index)=0;
     str=str(bool);
   end
end
%\********************************************************   

%/********************************************************
%linestyle
if isempty(str)
   stylecell=defaultstyle;
else
   stylecell=str;
end
%\********************************************************  

function[colors]=linestyle_colors
%COLORS
%
%  User-specified additional colors for use with LINESTYLE.
colors.k=[0 0 0];
colors.w=[1 1 1];
colors.b=[0 0 1];
colors.g=[0 0.5 0];
colors.r=[1 0 0];
colors.c=[0 0.75 0.75];
colors.m=[0.75 0 0.75];
colors.y=[0.75 0.75 0];
colors.A=(1-0)*[1 1 1];
colors.B=(1-1/10)*[1 1 1];
colors.C=(1-2/10)*[1 1 1];
colors.D=(1-3/10)*[1 1 1];
colors.E=(1-4/10)*[1 1 1];
colors.F=(1-5/10)*[1 1 1];
colors.G=(1-6/10)*[1 1 1];
colors.H=(1-7/10)*[1 1 1];
colors.I=(1-8/10)*[1 1 1];
colors.J=(1-9/10)*[1 1 1];
colors.K=(1-10/10)*[1 1 1];

colors.T=[     0    0.4470    0.7410];
colors.U=[0.8500    0.3250    0.0980];
colors.V=[0.9290    0.6940    0.1250];
colors.W=[0.4940    0.1840    0.5560];
colors.X=[0.4660    0.6740    0.1880];
colors.Y=[0.3010    0.7450    0.9330];
colors.Z=[0.6350    0.0780    0.1840];
    


