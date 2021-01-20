function[varargout]=printall(varargin)
%PRINTALL  Print and close all open figures.
%
%   PRINTALL prints all open figures to files in the current directory,
%   and closes the figures.  
%   
%   Files are named with the time stamp YYYY-MM-DD+HH|MM|SS followed
%   by their window number.
%
%   For image file output (png, jpeg, or tiff), the files are also cropped
%   using CROP by A. Bliss. 
%
%   PRINTALL(NAME) prepends the string NAME to the file names. 
%
%   PRINTALL(NAME,STR) passes the string STR as the printer device, for 
%   example, PRINTALL('-djpg -r200').  The default choice is STR='-dpng'.
%
%   PRINTALL(H,...) only prints those figures whose handles are given in H.
%
%   Usage: printall samplefigures
%          printall samplefigures -depsc
%          printall([],' -depsc')
%          printall(h,'-djpeg')
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2020 J.M. Lilly --- type 'help jlab_license' for details
 
str='-dpng';
name=[];
h=[];
if ~isempty(varargin)
    if ~ischar(varargin{1})
        h=varargin{1};
        varargin=varargin(2:end);
    end
end
if ~isempty(varargin)
    name=varargin{1};
    varargin=varargin(2:end);
end
if ~isempty(varargin)
    str=varargin{1};
    varargin=varargin(2:end);
end

if ~isempty(name)
    name=[name '@'];
end

ext=[];
if contains(str,'jpeg')
    ext='.jpg';
elseif contains(str,'png')
    ext='.png';
elseif contains(str,'tiff')
    ext='.tif';
end

if isempty(h)
    a=findobj('Type','Figure');
    for i=length(a):-1:1
        figure(a(i).Number)
        filename= datestr(now,31);filename(11)='+';filename(14:3:end)='|';
        namestr=[name filename '-' int2str(a(i).Number)];
        eval(['print ' str ' ' namestr])
        if ~isempty(ext)
            eval(['crop ' namestr ext])
        end
        close(a(i).Number);
    end
else
    for i=1:length(h)
        figure(h(i))
        filename= datestr(now,31);filename(11)='+';filename(14:3:end)='|';
        namestr=[name filename '-' int2str(h(i))];
        eval(['print ' str ' ' namestr])
        if ~isempty(ext)
            eval(['crop ' namestr ext])
        end
        close(h(i));
    end
end
%close all
