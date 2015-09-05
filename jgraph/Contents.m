% jGraph:  Plotting tools and refinements for making high-quality figures
%
% High-level graphical post-processing
%   fixlabels    - Specify precision of axes labels. 
%   fontsize     - Rapidly set title, axes, label, and text fontsizes.              
%   letterlabels - For automatically putting letter labels on subplots.           
%   linecolor    - Set line colors based on a property value within a colormap.                    
%   linestyle    - Rapidly set color, style, and width properties of lines.
%   linering     - Moves lines through the current line style order. 
%   packfig      - Squeeze together rows and/or columns of the current figure. 
%
%  Oceanography-specific plotting tools
%   hodograph  - Generate hodograph plots (simple and fancy).                                   
%   provec     - Generate progressive vector diagrams (simple and fancy).             
%   stickvect  - Plots "stick vectors" for multicomponent velocity time series.  
%   [See also CELLPLOT, REGIONPLOT, TOPOPLOT]
%
% Modified plotting functions
%   discretecolorbar - Plots a colorbar with discrete variation.  
%   jpcolor          - Modified version of PCOLOR appropriate for cell-centered grids.  
%   patchcontourf    - Generate filled contours using patches, with specified colors.
%   uvplot           - Plots the real and imaginary parts of a signal on the same axis.
%
% Figure tweaking
%   dlines     - Add diagonal lines to a plot.  
%   hlines     - Add horizontal lines to a plot.                                  
%   vlines     - Add vertical lines to a plot.   
%   nocontours - Removes contours from a CONTOURF plot.                           
%   noxlabels  - Remove some or all x-axis tick mark labels.                      
%   noylabels  - Remove some or all y-axis tick mark labels. 
%   xoffset    - Offsets lines in the x-direction after plotting.       
%   yoffset    - Offsets lines in the y-direction after plotting.                 
%   xtick      - Sets locations of x-axis tick marks.                             
%   ytick      - Sets locations of y-axis tick marks.                             
%   ztick      - Sets locations of z-axis tick marks.
%
% Graphics aliases
%   boxon      - Sets 'box' property to 'off'.                                    
%   boxon      - Sets 'box' property to 'on'.                                                                  
%   flipmap    - Flips the current colormap upside-down.                          
%   flipx      - Flips the direction of the x-axis                                
%   flipy      - Flips the direction of the y-axis    
%   inticks    - Sets the 'tickdir' property of the current axis to 'in'.
%   outticks   - Sets the 'tickdir' property of the current axis to 'out'.         
%   xlin       - Sets x-axis scale to linear.                                      
%   xlog       - Sets x-axis scale to logarithm.                                 
%   ylin       - Sets y-axis scale to linear.                                      
%   ylog       - Sets y-axis scale to log.        
%
% Colormaps
%  lansey     - The Lansey modification of Cynthia Brewer's "Spectral" colormap.
%
% Low-level functions
%   axeshandles - Returns handles to all axes children.
%   crop        - Gets rid of whitespace around an image. [by A. Bliss]             
%   linehandles - Finds all line and patch handles from a given set of axes.      
                
help jGraph
         
                                  
                             
