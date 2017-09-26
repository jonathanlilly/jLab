% jVarfun:  Perform common operations on multiple variables simultaneously.
%
%  Sizes and statistics
%   vsize      - Returns the sizes of multiple arguments.                           
%   vmean      - Mean over non-NaN elements along a specified dimension.             
%   vsum       - Sum over non-NaN elements along a specified dimension.              
%   vstd       - Standard deviation over non-NaN elements along a specfied dimension.
%   vmoment    - Central moment over non-NaN elements along a specfied dimension.    
%   vmedian    - Median over finite elements along a specified dimension. 
%
%  Reshaping, shifting, swapping
%   vcolon     - Condenses its arguments, like X(:).                                
%   vshift     - Cycles the elements of an array along a specified dimension.       
%   vsqueeze   - Squeezes multiple input arguments simultaneously.                  
%   vrep       - Replicates an array along a specified dimension.                   
%   vswap      - VSWAP(X,A,B) replaces A with B in numeric array X.                         
%   vtranspose - Transpose multiple input arguments simultaneously. 
%
%  Initializing
%   vzeros     - Initializes multiple variables to arrays of zeros or nans. 
%   vempty     - Initializes multiple variables to empty sets or empty cell arrays. 
%
%  Operations and indexing
%   vfilt      - Filtering along rows without change in length.                     
%   vdiff      - Length-preserving first central difference. 
%   vindex     - Indexes an N-D array along a specified dimension.                  
%   vindexinto - Indexes into N-D array along a specified dimension.
%
%  Reshaping dataset composed of irregular length segments 
%   col2mat    - Expands 'column-appended' data into a matrix.                      
%   colbreaks  - Insert NANs into discontinuties in a vector.                       
%   mat2col    - Compress NAN-padded matrix data into long columns.    
%   
%  See also jCell.
help jvarfun
               
          

