function[varargout]=jfig(varargin)
%JFIG  Shorthand for tweaking figures. 
%
%   JFIG is a shorthand syntax for tweaking figures.  It works like this.
%
%      jfig equal ticksout yrev axis|[5 45 5 45] caxis|[-8 8] ...
%           title|['Demonstration of JFIG'] colorbar|['eastoutside']
%
%   This call to JFIG is equivalent to
%     
%      axis equal, set(gca,'tickdir','out'), set(gca,'ydir','reverse')
%      axis([5 45 5 45]), caxis([-8 8]), title('Demonstration of JFIG')
%      colorbar(gca,'eastoutside')
%   
%   JFIG takes two kinds of arguments.  The first are one-word commands, 
%   such as 'equal', 'ticksout', and 'yrev'.  A list is provided below.  
%
%   The second kind of argument is a function with one input argument.  The
%   function name is separated from the argument by a vertical bar, '|'.
%   For these the function is simply called with the following argument.  
%
%   Note strings or arrays have to be enclosed in brackets in this format.
%   Arguments can also be variables, i.e., the following also works:
%
%       ax=[5 45 5 45]; cax=[-8 8];  jfig axis|ax caxis|cax
%
%   Any single-argument function can be called in this way.   
%
%   There are two exceptions for the single-argument functions.  Firstly, 
%   colorbar is called with GCA as its first argument whatever appears
%   after the '|' as its second argument.  Secondly, if the function name
%   is 'eval', its argument will be evaluated in the caller.  
%
%   The one-word commands and the single-argument functions can appear in
%   any order within a call to JFIG.
%
%   By default, JFIG turns on hold and also turns on the axis box. 
%   _______________________________________________________________________
%
%   One-word commands
%
%      new          - figure
%      clf          - clf
%      square       - axis square
%      equal        - axis equal
%      tight        - axis tight
%      holdon       - hold on
%      holdoff      - hold off
%      xlin         - set(gca,'xscale','linear')
%      ylin         - set(gca,'yscale','linear')
%      zlin         - set(gca,'zscale','linear')
%      xlog         - set(gca,'xscale','log')
%      ylog         - set(gca,'yscale','log')
%      zlog         - set(gca,'zscale','log')
%      xrev         - set(gca,'xdir','reverse')
%      yrev         - set(gca,'ydir','reverse')
%      zrev         - set(gca,'zdir','reverse')
%      xnorm        - set(gca,'xdir','normal')
%      ynorm        - set(gca,'ydir','normal')
%      znorm        - set(gca,'zdir','normal')
%      noxlabels    - set(gca,'xticklabel',[])
%      noylabels    - set(gca,'yticklabel',[])
%      nozlabels    - set(gca,'zticklabel',[])
%      xmonths      - set(gca,'xtick',[1:12]+0.5),xlabel('Month'),...
%                       set(gca,'xticklabel',['JFMAMJJASOND']')
%      ymonths      - set(gca,'ytick',[1:12]+0.5),ylabel('Month'),...
%                       set(gca,'yticklabel',['JFMAMJJASOND']')
%      zmonths      - set(gca,'ztick',[1:12]+0.5),zlabel('Month'),...
%                       set(gca,'zticklabel',['JFMAMJJASOND']')
%      ticksin      - set(gca,'tickdir','in')
%      ticksout     - set(gca,'tickdir','out')
%      boxon        - set(gca,'box','on')
%      boxoff       - set(gca,'box','off')
%      topoplot     - topoplot
%      nocontours   - nocnoturs
%      latlon       - xlabel('Longitude'),ylabel('Latitude')
%      landscape    - orient landscape
%      portrait     - orient portrait
%   _______________________________________________________________________
%
%   'jfig --f' generates a sample figure using the code presented above.
%
%   Usage: jfig equal ticksout yrev axis|[5 45 5 45] caxis|[-8 8] 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1},'--f')
   type makefigs_jfig
   makefigs_jfig;
   return
end

% str=[];
% for i=1:length(varargin)
%     %varargin{i}
%     str=[str ' ' varargin{i}];
% end

% varargin
% 
% for j=1:length(varargin)
%     str=varargin{j};
%     ii=min(strfind(str,' '));
%     if ~isempty(ii)
%         str=str(ii+1:end);
%     end
%     varargin{j}=str;
% end
% varargin

hold on, set(gca,'box','on')

for i=1:length(varargin)
    ii=findstr(varargin{i},'|');
    if ~isempty(ii)
        command=varargin{i}(1:ii-1);
        arg=varargin{i}(ii+1:end);
    else
        command=varargin{i};
        arg=[];
    end
    %command,arg
    if isempty(arg) 
        switch command
            %Expects zero
            case 'new',       figure
            case 'clf',       clf
            case 'square',    axis square
            case 'equal',     axis equal
            case 'tight',     axis tight
            case 'hold on',   hold on
            case 'hold off',  hold off
            case 'xlin',      set(gca,'xscale','linear')
            case 'xlog',      set(gca,'xscale','log')
            case 'ylin',      set(gca,'yscale','linear')
            case 'ylog',      set(gca,'yscale','log')
            case 'zlin',      set(gca,'zscale','linear')
            case 'zlog',      set(gca,'zscale','log')
            case 'xrev',      set(gca,'xdir','reverse')
            case 'yrev',      set(gca,'ydir','reverse')
            case 'zrev',      set(gca,'zdir','reverse')
            case 'xnorm',     set(gca,'xdir','normal')
            case 'ynorm',     set(gca,'ydir','normal')
            case 'znorm',     set(gca,'zdir','normal')
            case 'noxlabels', set(gca,'xticklabel',[])
            case 'noylabels', set(gca,'yticklabel',[])
            case 'nozlabels', set(gca,'zticklabel',[])
            case 'xmonths',   set(gca,'xtick',[1:12]+0.5),set(gca,'xticklabel',['JFMAMJJASOND']'),xlabel('Month')
            case 'ymonths',   set(gca,'ytick',[1:12]+0.5),set(gca,'yticklabel',['JFMAMJJASOND']'),ylabel('Month')
            case 'zmonths',   set(gca,'ztick',[1:12]+0.5),set(gca,'zticklabel',['JFMAMJJASOND']'),zlabel('Month')
            case 'ticksin',   set(gca,'tickdir','in')
            case 'ticksout',  set(gca,'tickdir','out')
            case 'boxon',     set(gca,'box','on')
            case 'boxoff',    set(gca,'box','off')
            case 'topoplot',  topoplot
            case 'nocontours',nocontours
            case 'latlon',    xlabel('Longitude'),ylabel('Latitude')
            case 'landscape', orient landscape
            case 'portrait',  orient portrait
        end
    else
        switch command
            case 'eval'
                evalin('caller',arg)
            case 'colorbar'
                evalin('caller',[command '(gca,' arg ');'])          
            otherwise
                %[command '(' arg ');']
                evalin('caller',[command '(' arg ');'])      
        end
    end
end


function[]=jfig_fragments
%These are old, unused fragments 
h=[];
packstr=[];
M=1;N=1;


%Get the first word
if strcmpi(str,'new')||strcmpi(str,'old')||strcmpi(str,'clf')
    %figure mode
    if strcmpi(str,'old')
        %if h is empty, look for subplots
        h=get(gcf,'children');
        bool=false(size(h));
        for i=1:length(h)
            if strcmpi(h(i).Type,'axes')
                bool(i)=true;
            end
        end
        h=h(bool);
    else
        if strcmpi(str,'new')
            figure
        elseif strcmpi(str,'clf')
            clf
        end
        %Subplot specification
        str=varargin{2};
        ii=strfind(str,'x');
        M=str2num(str(1:ii-1));
        N=str2num(str(ii+1:end));
        h=zeros(M*N,1);
        for i=1:M*N
            h(i)=subplot(M,N,i);
        end
        varargin=varargin(2:end);
    end
    varargin=varargin(2:end);
elseif strcmpi(str,'first')||strcmpi(str,'this')||strcmpi(str,'next')
    %axis mode, h is definitely a scalar
    if strcmpi(str,'first')
        h=getappdata(0,'jfigure_axeshandle');
        setappdata(0,'jfigure_currentaxes',1)
        axes(h(1));
        h=gca;
    elseif strcmpi(varargin{1},'next')
        h=getappdata(0,'jfigure_axeshandle');
        n=getappdata(0,'jfigure_currentaxes');
        n=n+1;
        if n>length(h)
            n=1;
        end
        setappdata(0,'jfigure_currentaxes',n)
        axes(h(n))
        h=gca;
    end
    varargin=varargin(2:end);
end

%At this point I have the axis handle I'll be working with 


%Look for figure-specific commands first
bool=true(size(varargin));
for j=1:length(varargin)
    str=varargin{j};
    if aresame(str(1:3),'pos')
        %Setting paper position
        ii=strfind(str,'|');
        str=str(ii+1:1:end);
        %str(real(str)==real('-'))=' ';
        eval(['set(gcf,''paperposition'',[' str '])'])
        bool(j)=false;
    elseif aresame(str(1:3),'pac')
        %Packing; save for later
        ii=strfind(str,'|');
        packstr=str(ii+1:1:end);
        bool(j)=false;
    elseif ~contains(str,'all')
        if contains(str,'x')&&(strfind(str,'x')~=1)&&(strfind(str,'x')~=length(str))
            %Subplot specification
            ii=strfind(str,'x');
            M=str2num(str(1:ii-1));
            N=str2num(str(ii+1:end));
            h=zeros(M*N,1);
            for i=1:M*N
                h(i)=subplot(M,N,i);
            end
            bool(j)=false;  
        end
    end
end
varargin=varargin(bool);
        
setappdata(0,'jfigure_axeshandle',h)
setappdata(0,'jfigure_currentaxes',length(h))

for i=1:length(h)
    axes(h(i)),hold on
    for j=1:length(varargin)
        str=varargin{j};
        if strcmpi(str(1),'x')||strcmpi(str(1),'y')||strcmpi(str(1),'z')
            if aresame(str(2),'|')
                %Limits are specified
                xyz=str(1);
                ii=strfind(str,'-');
                a=(str(3:ii-1));
                b=(str(ii+1:end));
                %['set(gca,' '''' str(1) 'lim'',[' a ',' b '])']
                eval(['set(gca,' '''' xyz 'lim'',[' a ',' b '])'])
            elseif aresame(str(2),'t')
                %Tickmarks are specified
                xyz=str(1);
                ii=strfind(str,'|');
                str=str(ii+1:1:end);
                str(real(str)==real('-'))=' ';
                eval(['set(gca,' '''' xyz 'tick'',[' str '])'])
            end
        else
            if aresame(str(1:3),'fli')
                %Flipping axis is specified
                xyz=str(5);
                %['set(gca,' '''' xyz 'dir'',''reverse'')']
                eval(['set(gca,' '''' xyz 'dir'',''reverse'')'])
            elseif aresame(str(1:3),'mon')
                %Monthly tick marks are specified
                xyz=str(end);
                str='[''JJFMAMJJASOND'']''';
                %['set(gca,' '''' xyz 'ticklabel'',[' str '])']
                eval(['set(gca,' '''' xyz 'tick'',[1.5:1:13.5])'])
                eval(['set(gca,' '''' xyz 'ticklabel'',[' str '])'])
            elseif aresame(str(1:3),'box')
                set(gca,'box','on')
            elseif aresame(str(1:3),'nob')
                set(gca,'box','off')
            elseif aresame(str(1:3),'all')
                %Apply to all subplots
                ii=strfind(str,'|');
                str=[str(ii(1)+1:1:end) ';'];
                evalin('caller',str)
%             elseif contains(str(1:3),'colorbar')
%                 %Add a colorbar
%                 2
%                 str
%                 ii=strfind(str,'|');
%                 str1=[str(ii(1)+1:1:end) ';'];
%                 hc=colorbar;
%                 hc.Label.String=str1;
            end
        end
    end
end
axes(h(1))
if ~isempty(packstr)
    packfig(M,N,packstr)
end

%After separating what is to the left and right of an option word
switch left
    case 'colorbar'
        h=color        
end

use dictionary
eval(command)

dictionary.xlog='set(gca,''xscale'',''log'')';
dictionary.xlin='set(gca,''xscale'',''lin'')';
dictionary.ylog='set(gca,''yscale'',''log'')';
dictionary.ylin='set(gca,''yscale'',''lin'')';
dictionary.zlog='set(gca,''zscale'',''log'')';
dictionary.zlin='set(gca,''zscale'',''lin'')';
dictionary.zlin='colormap(flipud(colormap))';
dictionary.inticks='set(gca,''tickdir'',''in'')';
dictionary.outticks='set(gca,''tickdir'',''out'')';
dictionary.box='set(gca,''box'',''on'')';
dictionary.nobox='set(gca,''box'',''off'')';
%dictionary.xflip=['if strcmpi(get(gca,''xdir''),''normal'')' setstr(10)...
%       set(gca,'xdir','reverse')  setstr(10)...
% elseif strcmpi(get(gca,'xdir'),'reverse') 
%   set(gca,'xdir','normal')
% end

command='hey';

        
dictionary.outticks='set(gca,''tickdir'',''out'')';
dictionary.box='set(gca,''box'',''on'')';
dictionary.nobox='set(gca,''box'',''off'')';


%function[]=jfig_test
%reporttest('JFIGURE',aresame())
%jfig new 3x2 x|4-10 xtick|[1.5:1:12.5] flipy monthsx box
