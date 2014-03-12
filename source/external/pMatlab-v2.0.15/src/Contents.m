% pMatlab: Parallel Matlab Toolbox v2.0
%
%   For more information on any of the functions listed below, type 
%       help <function_name> 
%   at the MATLAB command prompt. (Note: please include 'dmat/' or 
%   'map/' before the function name to get help for the
%   appropriate overloaded function.)
%   
%   ---------- User Functions/Scripts ----------
%   Np                      Number of Matlab instances launched.
%   Pid                     ID of current Matlab instance.
%   pRUN                    Command for launching a pMatlab script.
%   pMatlab_ver             Returns the pMatlab version as a string.
%   RecvMsg                 Receive a message from another Pid.
%   SendMsg                 Send a message to another Pid.
%   synch                   No-op if X is a DOUBLE.
%   dmat/abs                Absolute value.   
%   dmat/acos               Inverse cosine.
%   dmat/acosd              Inverse cosine, result in degrees.
%   dmat/acosh              Inverse hyperbolic cosine.
%   dmat/acot               Inverse cotangent.
%   dmat/acotd              Inverse cotangent, result in degrees.
%   dmat/acoth              Inverse hyperbolic cotangent.
%   dmat/acsc               Inverse cosecant.
%   dmat/acscd              Inverse cosecant, result in degrees.
%   dmat/acsch              Inverse hyperbolic cosecant.
%   dmat/angle              Phase angle.
%   dmat/asec               Inverse secant.
%   dmat/asecd              Inverse secant, result in degrees.
%   dmat/asech              Inverse hyperbolic secant.
%   dmat/asin               Inverse sine.
%   dmat/asind              Inverse sine, result in degrees.
%   dmat/asinh              Inverse hyperbolic sine.
%   dmat/atan               Inverse tangent.
%   dmat/atand              Inverse tangent, result in degrees.
%   dmat/atanh              Inverse hyperbolic tangent.
%   dmat/ceil               Round towards plus infinity.
%   dmat/complex            Construct complex distributed matrix from real 
%                           distributed matrix.
%   dmat/conj               Complex conjugate.
%   dmat/cos                Cosine.
%   dmat/cosd               Cosine of argument in degrees.
%   dmat/cosh               Hyperbolic cosine.
%   dmat/cot                Cotangent.
%   dmat/cotd               Cotangent of argument in degrees.
%   dmat/coth               Hyperbolic cotangent.
%   dmat/csc                Cosecant.
%   dmat/cscd               Cosecant of argument in degrees.
%   dmat/csch               Hyperbolic cosecant.
%   dmat/display            Display distributed array.            
%   dmat/double             Convert each local part of the DMAT to double
%                           precision.
%   dmat/eq                 == Equal. 
%   dmat/exp                Exponential.
%   dmat/expm1              Compute exp(x)-1 accurately.
%   dmat/fix                Round towards zero.
%   dmat/floor              Round towards minus infinity.
%   dmat/imag               Complex imaginary part.
%   dmat/int16              Convert each local part of the DMAT to signed
%                           16-bit integers.
%   dmat/int32              Convert each local part of the DMAT to signed
%                           32-bit integers.
%   dmat/int64              Convert each local part of the DMAT to signed
%                           64-bit integers.
%   dmat/int8               Convert each local part of the DMAT to signed
%                           8-bit integers.
%   dmat/log                Natural logarithm.
%   dmat/log1p              Compute log(1+x) accurately.
%   dmat/log10              Common (base 10) logarithm.
%   dmat/log2               Base 2 logarithm.
%   dmat/pow2               Base 2 power.
%   dmat/real               Complex real part.
%   dmat/realpow            Power that will error out on complex result.
%   dmat/reallog            Natural logarithm of real number.
%   dmat/realsqrt           Square root of number greater than or equal to
%                           zero.
%   dmat/round              Round towards nearest integer.
%   dmat/sec                Secant.
%   dmat/secd               Secant of argument in degrees.
%   dmat/sech               Hyperbolic secant.
%   dmat/sign               Signum.
%   dmat/sin                Sine.
%   dmat/sind               Sine of argument in degrees.
%   dmat/single             Convert each local part of the DMAT to single
%                           precision.
%   dmat/sinh               Hyperbolic sine.
%   dmat/size               Size of the distributed array.
%   dmat/sparse             Converts a distributed matrix to a sparse
%                           distributed matrix
%   dmat/sqrt               Square root.
%   dmat/subsasgn           Subscripted assignment to a distributed object.
%                           Called for syntax A(I) = B.
%                           Should not be called directly.
%   dmat/subsref            Subscripted reference. Called for syntax A(I) = B.
%   dmat/synch              Syncronize the data in the distribute matrix. 
%   dmat/tan                Tangent.
%   dmat/tand               Tangent of argument in degrees.
%   dmat/tanh               Hyperbolic tangent.
%   dmat/uint16             Convert each local part of the DMAT to unsigned
%                           16-bit integers.
%   dmat/uint32             Convert each local part of the DMAT to unsigned
%                           32-bit integers.
%   dmat/uint64             Convert each local part of the DMAT to unsigned
%                           64-bit integers.
%   dmat/uint8              Convert each local part of the DMAT to unsigned
%                           8-bit integers.
%   map/display             Display map object.
%   map/eq                  == Equal.
%   map/map                 Map class constructor.             
%   map/ne                  ~= Not equal.
%   map/ones                Ones distributed array.
%   map/rand                Rand distributed array.
%   map/spalloc             Sparse distributed array.
%   map/subsasgn            Subscripted assignment. Should not be called
%                           directly.
%   map/subsref             Subscripted reference. Should not be called
%                           directly.
%   map/zeros               Zeros distributed array.
%
%
%   ---Please use the following with caution, since implementation is only
%       provided for some special cases.
%   transpose_grid          Redistributes a dmat by transposing its grid.
%   dmat/conv2              Two dimensional convolution.
%   dmat/ctranspose         ' Complex conjugate transpose.
%   dmat/dct                Distributed discrete cosine transform.
%   dmat/fft                Discrete Fourier transform on a distributed matrix.
%   dmat/find               Find indices of nonzero elements.
%   dmat/gt                  > Greater than.
%   dmat/ldivide            .\ Left array divide.
%   dmat/minus              - Minus.
%   dmat/mtimes             * Matrix multiply.
%   dmat/plus               + Plus.
%   dmat/power              .^ Array power.
%   dmat/rdivide            ./ Right array divide.
%   dmat/summation          Takes a sum of several DMATs.
%   dmat/times              .* Array multiply.
%   dmat/transpose          .' Transpose.
%
%
%   ---------- Library Functions That a User Might Use ----------
%   binfun                  Binary operation for non-distributed arrays.
%   global_block_range      Returns the range of indices in the specified
%                           dimension.
%   global_block_ranges     Returns the range of indices in the specified
%                           dimension.
%   global_range            Returns the range of indices in the specified
%                           dimension.
%   global_ranges           Returns the range of indices in the specified
%                           dimension.
%   global_ind              Returns the indices in the specified dimension.
%   global_inds             Returns indices in the specified dimension.
%   g_mat                   No-op if the input argument is a DOUBLE. 
%   local                   No-op if the input argument is a DOUBLE.
%   put_local               Mimics dmat/put_local on a DOUBLE.
%   remap                   Remaps a distributed array.
%   dmat/binfun             Binary elementwise operation.
%   dmat/global_block_range Returns the ranges of global indices local to the
%                           current processor.
%   dmat/global_block_ranges Returns the ranges of global indices for all
%                           processors in the map of distributed array D on 
%                           all processors in communication scope.
%   dmat/global_ind         Returns the global indices local to the
%                           current processor. 
%   dmat/global_inds        Returns the global indices for all
%                           processors in the map of distributed array D.
%   dmat/global_range       Returns the ranges of global indices local to the
%                           current processor.
%   dmat/global_ranges      Returns the ranges of global indices for all
%                           processors in the map of distributed array D.
%   dmat/local              Returns the local part of the distributed array.
%   dmat/ndims              Returns the number of dimensions of the DMAT.
%   dmat/put_local          Assigns new data to the local part of the
%                           distributed array.
%   map/sparse              Sparse distributed array.  SPALLOC is the
%                           recommended method of creating sparse
%                           distributed arrays.
%
%
%   ---------- Library Functions ----------
%   agg                     No-op if the input argument is a DOUBLE.
%   agg_all                 No-op if the input argument is a DOUBLE.
%   exists_falls_intersection Returns TRUE is the intersection of two FALLS is
%                           non-empty, FALSE otherwise.
%   falls_intersection      Given two FALLS, find the intersection FALLS.
%   gen_falls               Generates the local falls structure for processor P
%                           from the given PITFALLS structure.
%   gen_pitfalls            Given the number of processors, distribution spec,
%                           and the length of the dimension, generates all the
%                           PITFALLS information.
%   get_global_ind          Returns a cell array of global indices stored
%                           locally given a FALLS object 
%   get_ind_range           Creates index ranges.
%   get_local_falls         Given the PITFALLS object, the grid and local
%                           processor rank, returns an array of local FALLS
%                           objects.
%   get_local_ind           Returns a cell array of local indices given an array
%                           of global indices on the current processor and an
%                           array of global indices that are being referenced.
%   get_local_proc          Returns the rank of the processor that contains
%                           index IND.
%   local_dims              Given the local FALLS and object dimension,
%                           calculates a vector of local dimensions of the
%                           distributed object.
%   ls_intersection         Finds an intersection of two line segments.
%   multicast               Sends data from a single source to multiple
%                           destinations.
%   n_dim_find              N-dimensional FIND functions.
%   pMatlabGlobalsInit      Initializes global variables.
%   pRUN_Parallel_Wrapper   Wrapper script for launching pMatlab scripts.
%   pMatlab_Init            *Deprecated* Initializes all the needed MPI variables.
%   pMatlab_Finalize        *Deprecated* Exits non-leader MATLAB processes.
%   pitfalls_intersection   Computes intersection of two PITFALLS.
%   reconstruct             Given collected distributed data in grid layout,
%                           reconstructs the original object according to data
%                           distribution.
%   resize                  Helper function used by distributed mtimes.
%   dmat/agg                Aggregates the parts of a distributed matrix on the
%                           leader processor.
%   dmat/agg_all            Aggregates the parts of a distributed matrix on all
%                           processors in the map of the distributed matrix
%   dmat/dmat               Distributed matrix constructor.
%   dmat/g_mat              Aggregate parts of the distributed matrix.
%   dmat/grid               Returns the processor grid onto which the
%                           distributed array is mapped.
%   dmat/pitfalls           Returns the pitfalls structure associated with the
%                           distributed array.
%   dmat/submat             Helper function used by subsref.
%   dmat/subsasgn_data      Helper function for distributed array subsasgn. 
%   dmat/subsasgn1D         One dimensional subsasgn.
%   dmat/subsasgn2D         Two dimensional subsasgn.
%   dmat/subsasgn3D         Three dimensional subsasgn.
%   dmat/subsasgn4D         Four dimensional subsasgn.
%   map/inmap               Checks if a processor is in the map.
%   map/transpose           Transposes a map.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMatlab: Parallel Matlab Toolbox
% Software Engineer: Ms. Nadya Travinin (nt@ll.mit.edu)
% Architect:      Dr. Jeremy Kepner (kepner@ll.mit.edu)
% MIT Lincoln Laboratory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2005, Massachusetts Institute of Technology All rights 
% reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are 
% met:
%      * Redistributions of source code must retain the above copyright 
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright 
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the
%        distribution.
%      * Neither the name of the Massachusetts Institute of Technology nor 
%        the names of its contributors may be used to endorse or promote 
%        products derived from this software without specific prior written 
%        permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
