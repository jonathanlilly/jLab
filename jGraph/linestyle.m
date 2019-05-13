function[]=linestyle(varargin)
%LINESTYLE  Rapidly set color, style, and width properties of lines.
%	    	    	    
%   LINESTYLE provides an efficient way to set common properties of large
%   groups of lines, and allows the user to quickly apply predefined sets 
%   of line colors, styles, and widths designed for different purposes.
%
%   'LINESTYLE 1k-- 2.5g- m-.' sets the first line to width 1, color black,
%   and dashed style; the second to width 2.5, color green, and solid 
%   style; and the third to width 1, color magenta, and dash-dotted style.
%
%   LINESTYLE's color options supports grayscale colors, denoted by capital
%   letters A--F, with A being white and K being black.  
%
%   Also, the new Matlab colors are denoted by the letters T through Z.
%
%   LINESTYLE supports a number of predefined formats. Type LINESTYLE with
%   no arguments to see a list of current style sets. 
%  
%   In general the format is 
%
%       LINESTYLE STR1 STR2 STR3 ...   
%
%   where each of the STRs may contain a number, specifying the width of 
%   the Ith line; a letter, for the color; and a style string.  Any two of
%   these are optional, with a unit width solid black line as the default.
%
%   By default, styles are looped if the number of STRs input is less than
%   the number of lines in the current plot. The input
%
%       LINESTYLE STR1 STR2 ... STRN +++
%
%   causes the last style input, STRN, to be repeated instead.      
%
%   Note that LINESTYLE will set styles for both lines as well as patch
%   edges.  This is usually the desired behavior, as contours from CONTOUR
%   are actually patches.  If this is not what is desired, call LINESTYLE
%   with the explicit handle specification below, or else add contour plots
%   to the current plot after calling LINESTYLE.
%   _______________________________________________________________________
%
%   Locking and unlocking 
%
%   LINESYTLE LOCK and LINESTYLE UNLOCK lock and unlock all axes in the 
%   current figure.  When LOCK is on, calls to LINESTYLE or LINERING are 
%   applied to all lines in the current figure.
%   _______________________________________________________________________
%
%   Handle specification
%  
%   LINESTYLE -H HAN STR1 STR2 .... applies the formatting only to the line
%   handles contained in handle array HAN.
%
%   LINESTYLE(HAN,STR) also works.
%
%   HAN may also be the handle to a group of contours.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2016 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==1
    if strcmpi(varargin{1},'--t')
        return
    end
end

h=[];
linestyles=linestyle_linestyles;
bextendlast=0;
bloop=1;
bhandleinput=0;
na=nargin;
num=1;

nflagsin=0;
hflagin=0;

vars=varargin;

if na==0
   %if no arguments, just display available linestyle sets
   disp(' ');
   disp('Available linestyle sets:');
   disp(' ');
   disp(linestyles)
   disp(' ');
   return
end


if na==1
   %Special cases 
   if strcmpi(vars{1},'lock')
        setappdata(gcf,'linestylelock',1)
	return
   elseif strcmpi(vars{1},'unlock')
        setappdata(gcf,'linestylelock',0)
	return
   end
end

if ischar(vars{1})
    if length(vars{1})>1
        if strcmpi(vars{1}(1:2),'-h')
           %-h flag for handle input
           hflagin=1;
           nflagsin=1;

           bhandleinput=1;
           h=vars{2};
           if ischar(h)
               eval(to_grab_from_caller(2))
               eval(['h=' varargin{2} ';']);
           end
           %h=flipud(h(:));
           na=na-2;       
           vars=vars(3:end);
        end
    end
else
    if na==2
        h=vars{1};
        bhandleinput=1;
        na=na-1;
        vars=vars(2:end);
    end
end
        
axes_original=gca;

if isempty(h)&&bhandleinput 
    return;  %Do nothing; we've been given an empty array of handles
end
if isempty(h)
    h=flipud(linehandles(gca));
else
    %Change to axes to which the given handles belong
    hparent=get(h(1),'parent');
    if strcmpi(get(hparent,'type'),'axes')
        axes(hparent)
    else
        axes(get(hparent,'parent'))
    end
end

if length(h)==1
    if strcmpi(get(h,'type'),'hggroup')
        h=get(h,'children');
    end
end

if na==1
    if isfield(linestyles,vars{1})
        [widthcell,stylecell,colorcell]=linestyleapply(h,linestyles,vars{1});
        bloop=0;  %supress loop if name input
    end
end

 
% %/********************************************************
% %Block for applying LINESTYLE to input line handle
% if nargin==2
%   if ~ischar(vars{1});
%     h=vars{1};
%     sty=vars{2};
%     if isfield(linestyles,sty)
%        [widthcell,stylecell,colorcell]=linestyleapply(linestyles,sty);
%        bloop=0;  %supress loop if name input
%     else
%        [widthcelli,stylecelli,colorcelli]=linestyleparse(sty);    
%        widthcell{1}=widthcelli;
%        stylecell{1}=stylecelli;
%        colorcell{1}=colorcelli;
%     end
%     
%     bextendlast=1;
%     bhandleinput=1;
%     %bloop=0;  %supress loop if line handles are input
%   end
% end
% %\********************************************************
%    
% 



if isappdata(gca,'lineringpointer')
  num=getappdata(gca,'lineringpointer');   %Remember current linering
end

if bloop
   %This is specifically for when called like linestyle(h,'2g 1g 1r')
   %since normally Matlab would see this as one string, not three
   if na==1
      temp=vars{1};
      indexspaces=strfind(temp,' ');
      if ~isempty(indexspaces)
         a=[1 indexspaces+1];
         b=[indexspaces-1 length(temp)];
         vars=[];
         jj=0;
         for i=1:length(a)
             if ~isempty(temp(a(i):b(i)))
                 jj=jj+1;
                vars{i}=temp(a(i):b(i));
             end
         end
         na=length(vars);
      end
   end
   for i=1:na
     temp=vars{i};
     if strcmpi(temp,'+++')  %ending in continuation symbol
         bextendlast=1;
     else
        [widthcelli,stylecelli,colorcelli]=linestyleparse(temp);
        widthcell{i}=widthcelli;
        stylecell{i}=stylecelli;
        colorcell{i}=colorcelli;
     end  
   end
end

N=length(h);
M=length(stylecell);
if bextendlast
  %Extend by looping continuing last
  for j=M+1:N
      colorcell{j}=colorcell{M}; 
      stylecell{j}=stylecell{M}; 
      widthcell{j}=widthcell{M}; 
  end
else
  %Extend by looping whole structure 
  while N>M
      colorcell((1:M)+M)=colorcell;
      stylecell((1:M)+M)=stylecell; 
      widthcell((1:M)+M)=widthcell;  
      M=length(stylecell);
  end
end

%Truncate if too long
colorcell=colorcell(1:N)';
stylecell=stylecell(1:N)';
widthcell=widthcell(1:N)';

locked=0;      %Check to see if current figures axes are locked
if isappdata(gcf,'linestylelock')
   if getappdata(gcf,'linestylelock')
        locked=1;	
   end
end

[boolline,boolpatch]=handletype(h);

if ~locked || bhandleinput    %Just apply to current axes
  linering(1*sqrt(-1));%Go to first position in linering
  set(h(boolline),{'color'},colorcell(boolline),{'linewidth'},widthcell(boolline),{'linestyle'},stylecell(boolline))
  set(h(boolpatch),{'edgecolor'},colorcell(boolpatch),{'linewidth'},widthcell(boolpatch),{'linestyle'},stylecell(boolpatch))
  linering(num*sqrt(-1));%Return to original
elseif locked && (~bhandleinput)
  h1=axeshandles(gcf);	
  for i=1:length(h1)   %Loop over all axes
       linering(1*sqrt(-1));%Go to first position in linering
       h=linehandles(h1(i));
       [boolline,boolpatch]=handletype(h);
       set(h(boolline),{'color'},colorcell(boolline),{'linewidth'},widthcell(boolline),{'linestyle'},stylecell(boolline))
       set(h(boolpatch),{'edgecolor'},colorcell(boolpatch),{'linewidth'},widthcell(boolpatch),{'linestyle'},stylecell(boolpatch))
       linering(num*sqrt(-1));%Return to original
  end
end

axes(axes_original)

function[boolline,boolpatch]=handletype(h)

htype=get(h,'type');
if ~iscell(htype)
    htypetemp=htype;
    clear htype
    htype{1}=htypetemp;
end

boolline=false(size(h));
boolpatch=false(size(h));
for j=1:length(h)
    boolline(j)= strcmpi(htype{j},'line');
    boolpatch(j)= strcmpi(htype{j},'patch');   
end


function[widthcell,stylecell,colorcell]=linestyleapply(h,linestyles,name)
%Apply contents of LINESTYLE definitions to cell arrays
  
% %account for the fact that LINERING flips the linering
% %to put the first line on the top
% if isappdata(gca,'lineringflipped')
%    if getappdata(gca,'lineringflipped')
% 	h=flipud(h);
%    end
% end
  
colors=linestyleparse;
%h=flipud(linehandles(gca));

colorx=[];
stylex=[];
widthx=[];


if ~isfield(linestyles,name)
  error('That is not a valid linestyle set name.')
else
  colorx=getfield(linestyles,name,{1});
  colorx=colorx{1};
  stylex=getfield(linestyles,name,{2});
  stylex=stylex{1};
  widthx=getfield(linestyles,name,{3});
  widthx=widthx{1};
end

%make the vectors the same length as the handles by cycling
brepeat=0;
N=length(colorx)-3;
if ~isempty(colorx)
   colorx=colorx(:);
   if length(colorx)>3
	if strcmpi(colorx(end-2:end)','...')
	     colorx(end-2:length(h),:)=colorx(end-3,:);
	     brepeat=1;
	end
   end 		    
   while size(colorx,1)<length(h)
	 colorx=[colorx;colorx];
   end
   colorx=colorx(1:length(h),:);
end

if ~isempty(stylex)
   if size(stylex,1)==1
      stylex=stylex(:);
   end
   if length(stylex)==N && brepeat
      stylex(end:length(h),:)=stylex(end,:);
   else
       while size(stylex,1)<length(h)	
           stylex=[stylex;stylex];
       end
   end
   stylex=stylex(1:length(h),:);
end

if ~isempty(widthx)
   widthx=widthx(:);
   if length(widthx)==N && brepeat
      widthx(end:length(h),:)=widthx(end,:);
   else
       while size(widthx,1)<length(h)	
           widthx=[widthx;widthx];
       end
   end
   widthx=widthx(1:length(h),:);
end

%convert color into numeric values
temp=colorx;
clear colorx
colorx=[];
for i=1:length(temp)
    colorx(i,:)=getfield(colors,temp(i));
end


%now put into cell arrays, because that's what matlab understands
colorcell=[];
stylecell=[];
widthcell=[];

for i=1:length(h)
    if ~isempty(colorx)
	colorcell{i}=colorx(i,:);
    end    
    if ~isempty(stylex)
	stylecell{i}=stylex(i,:);
    end
    if ~isempty(widthx)
	widthcell{i}=widthx(i);
    end
end

%kludge
%h=flipud(h);
if ~isempty(widthcell),set(h,{'linewidth'},widthcell'),end
if ~isempty(stylecell),set(h,{'linestyle'},stylecell'),end
if ~isempty(colorcell),set(h,{'color'},colorcell'),end


function[linestyles]=linestyle_linestyles
%LINESTYLES
%
%  User-specified line style sets for use with LINESTYLE.
%
%    linestyles.name={'COLOR','STYLE','WIDTH'};
%
%  Line handles will be cycled through the available styles,
%  colors, and widths, so that S, C, and W need not have the 
%  same length as the number of line handles nor as each other.
%
%  An ellipsis may be used (for C only) to indicate a repeated
%  last color, thus 'kbgr...' is equivalent to 'krbrrrrrrr...'.

    
linestyles.colors={'bgrcmy','-',1};
linestyles.colors2={'bgrcmy','-',2};
linestyles.colors3={'bgrcmy','-',3};
linestyles.colors4={'bgrcmy','-',4};
%linestyles.default={'bgrcmyk','-',1};
%linestyles.default2={'bgrcmyk','-',2};
%linestyles.default3={'bgrcmyk','-',3};
%linestyles.default4={'bgrcmyk','-',4};
linestyles.former={'bgrcmyk','-',1};
linestyles.default={'TUVWXYZ','-',1};
linestyles.default2={'TUVWXYZ','-',2};
linestyles.default3={'TUVWXYZ','-',3};
linestyles.default4={'TUVWXYZ','-',4};
linestyles.black={'k','-',1};
linestyles.just1={'bDDDDDD...','-',[2 1 1 1 1 1 1]};
linestyles.just2={'bgDDDDD...','-',[2 2 1 1 1 1 1]};
linestyles.just3={'bgrDDDD...','-',[2 2 2 1 1 1 1]};
linestyles.just4={'bgrcDDD...','-',[2 2 2 2 1 1 1]};
linestyles.just5={'bgrcmDD...','-',[2 2 2 2 2 1 1]};
linestyles.just6={'bgrcmyD...','-',[2 2 2 2 2 2 1]};
linestyles.groups4={'kkkkbbbbggggrrrr','-',1};
linestyles.groups5={'kkkkkbbbbbgggggrrrrr','-',1};
linestyles.groups6={'kkkkkkbbbbbbggggggrrrrrr','-',1};
linestyles.groups7={'kkkkkkkbbbbbbbgggggggrrrrrrr','-',1};
linestyles.groups8={'kkkkkkkkbbbbbbbbggggggggrrrrrrrr','-',1};
linestyles.groups9={'kkkkkkkkkbbbbbbbbbgggggggggrrrrrrrrr','-',1};
linestyles.groups10={'kkkkkkkkkkbbbbbbbbbbggggggggggrrrrrrrrrr','-',1};
linestyles.thick={'TUVWXYZ','-',2};
linestyles.two={'Ek',['- ';'--'],[2 1]};
linestyles.three={'kEk',['--';'- ';'- '],[2 2 1]};
linestyles.four={'EkEk',['--';'--';'- ';'- '],[2 1 2 1]};
linestyles.just2bw={'kkkkkkk...','-',[2.5 2.5 1 1 1 1 1]};
linestyles.triplets={'bbbgggrrr',['- ';'--';'- '],[1 2 2]};
linestyles.redgreenblue={'rgb',['-';'-';'-'],[1 1 1]};
linestyles.dots={'bgrcmyk','.',1};
linestyles.graydots={'D','.',1};

