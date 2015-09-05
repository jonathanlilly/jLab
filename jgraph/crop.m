function crop(filename,append,margin)
%CROP  Gets rid of whitespace around an image. [by A. Bliss]
%
%   CROP(FILENAME,APPEND,MARGIN) is the full calling form. APPEND and
%      MARGIN are optional inputs.
%
%   CROP('filename.ext') crops the image in the file and saves it using
%   the original filename, overwriting the old image. The extension (ext)
%   can be anything IMREAD supports.
%
%   CROP(directory) crops all images in a directory. 
%
%   If APPEND is 1, CROP saves the cropped image as 'filename_crop.ext'
%   in the same directory as the original.
%
%   MARGIN sets the margin width in pixels (default is 10).
%
%   Changes since version 1:
%       1. Now accepts directories as input.
%       2. Improved input checks and comments.
%       3. Handles transparency in .png files.
%       4. Just prior to saving version 1, I added the APPEND option.
%
%   Example:
%       crop('C:\MATLAB7\toolbox\matlab\demos\html\cruller_01.png',1)
%
%   See also: IMREAD, IMWRITE
%   See crop_license.m for the license details.
%   Written by Andy Bliss Sept 8th, 2006. Revised May 31, 2012.
%   __________________________________________________________________
%   This third-party software is included with JLAB for convenience, in
%   accordance with the redistribution policy at crop_license.m.

%   Requirements: your FIND must allow the 'last' option (version 7+?)

%set the default margin width in pixels
if nargin<3 || isempty(margin)
    margin=10; %15;
end
%default is to save with original filename
if nargin<2 || isempty(append)
    append=0;
end

%Get image names
if isstruct(filename) %assuming it is a struct as produced from DIR
    files=filename;
elseif isdir(filename) %if the input is a directory, get all the image files from the directory
    currentdir=pwd;
    cd (filename)
    files=[dir('*.png'); dir('*.gif'); dir('*.bmp'); dir('*.jpg')];
else %if it is a single file:
    files.name=filename;
end

%loop over all the files
for n=1:length(files)
    filename=files(n).name;

    %get file info
    info=imfinfo(filename);
    
    %get the image
    if strcmpi(info.Format,'png') %if it has transparent pixels
        T=imread(filename,'backgroundcolor',[1 1 1]); %backgroundcolor makes transparent pixels white, so they don't affect cropping.
    else
        T=imread(filename);
    end

    %sum the RGB values of the image
    xsum=sum(sum(T,3));
    ysum=sum(sum(T,3),2);
    % figure,plot(xsum),title('xsum'),xlabel('distance from left (pixels)'),ylabel('image intensity (big numbers are white)')
    % figure,plot(ysum),title('ysum'),xlabel('distance from top (pixels)'),ylabel('image intensity (big numbers are white)')

    %xsum will be equal to max(xsum) wherever there is a blank column in 
    %   the image (rgb white is [255,255,255]). The left edge for the 
    %   cropped image is found by looking for the first column in which 
    %   xsum is less than max(xsum) and then subtracting the margin.
    %   Similar code for other edges.
    xleftedge=find(xsum<max(xsum),1,'first')-margin;
    if xleftedge<1
        xleftedge=1;
    end
    xrightedge=find(xsum<max(xsum),1,'last')+margin;
    if xrightedge>length(xsum)
        xrightedge=length(xsum);
    end
    ytopedge=find(ysum<max(ysum),1,'first')-margin;
    if ytopedge<1
        ytopedge=1;
    end
    ybottomedge=find(ysum<max(ysum),1,'last')+margin;
    if ybottomedge>length(ysum)
        ybottomedge=length(ysum);
    end

    %resave the image
    if append
        filename=[filename(1:end-4) '_crop' filename(end-3:end)];
    end
    imwrite(T(ytopedge:ybottomedge,xleftedge:xrightedge,:),filename)
end
%change back to calling directory, if necessary
if exist('currentdir','var')
    cd(currentdir)
end
