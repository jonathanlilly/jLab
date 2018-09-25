% jCell:  Tools for operating on cell arrays of column vectors
%
% Basic mathematical operations
%   cellabs    - Absolute value of each element in a cell array.                    
%   cellmax    - Maximum of each element in a cell array.                           
%   cellmin    - Minimum of each element in a cell array. 
%   cellsum    - Sum of each element a cell array, possibly weighted.
%   cellmean   - Mean value of each element a cell array, possibly weighted.  
%   cellstd    - Standard deviation of each element a cell array, possibly weighted.
%   cellmed    - Median value of each element a cell array.
%   cellreal   - Real part of each element in a cell array.                         
%   cellimag   - Imaginary part of each element in a cell array.    
%   celllog10  - Base ten logarithm of each element in a cell array.
%   celladd    - Addition acting on each element in a cell array.                   
%   cellmult   - Multiplication acting on each element in a cell array.   
%   celldiv    - Division acting on each element in a cell array.
%
% Reshaping, indexing, and sizes
%   cell2col   - Converts cell arrays of numeric arrays into 'column-appended' form.
%   col2cell   - Converts 'column-appended' data into cell arrays of numeric arrays.
%   cellindex  - Applies a cell array of indices to a cell array of column vectors. 
%   cellchunk  - Converts cell array data into uniform length 'chunks'.             
%   cellength  - Length of each element in a cell array.                            
%   cellsize   - Size of each element in a cell array along specified dimension.
%   cellget    - Indexes a cell array of numerical arrays by ID number.
%   cellimit   - Limits the ranges of times in a cell array of numerical arrays.
%
% Data processing
%   cellstrip  - Strips INF values from the beginnings or ends of cell arrays.
%   cellsplit  - Splits cell arrays of numeric arrays at data gaps.
%   cellprune  - Removes all empty cells, or cells less than a specified length.
%   cellpack   - Removes all INF values from cell arrays of numeric arrays.
%   cellfill   - Fills missing data marked by INFs in a cell array.
%   cellgrid   - Interpolate a cell array of numeric arrays onto a regular grid.
%   cellfirst  - Returns the first element of each entry in a cell array.
%
% Plotting
%   cellplot   - Rapidly plot all element of a cell array of numeric arrays. 
%
% See also jVarfun.
help jCell
