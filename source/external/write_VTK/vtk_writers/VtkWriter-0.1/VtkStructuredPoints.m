% Copyright (c) 2012, Michael Mallon
% All rights reserved.

% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:

% Redistributions of source code must retain the above copyright notice, this list
% of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, this
% list of conditions and the following disclaimer in the documentation and/or
% other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

classdef VtkStructuredPoints < handle
    properties (SetAccess=private)
        n
        org
        spacing
    end
    methods
        function obj = VtkStructuredPoints(varargin)
            switch length(varargin)
                case 0
                case 3
                    if (all(size(varargin{1})==size(varargin{2})) && ...
                            all(size(varargin{1})==size(varargin{3})) && ...
                            (all(size(varargin{1})==[1 1]) || ...
                             all(size(varargin{1})==[1 2]) || ...
                             all(size(varargin{1})==[1 3]) ))
                        obj.n = varargin{1};
                        obj.org = varargin{2};
                        obj.spacing = varargin{3};
                    else
                        throw(MException('MATLAB:argchk',...
                        'All parameters must have the same size and/or must be 1D, 2D or 3D.'));
                    end
                otherwise
                    throw(MException('MATLAB:nargchk',...
                        ['Structured Points constructor takes 0 or 3 arguments:\n'...
                        'VtkStructuredPoints() or VtkStructuredPoints(n,origin,spacing)']));
            end
        end
        
        function obj = setNumPoints(obj,n)
            obj.n = n;
        end
        
        function obj = setOrigin(obj,orgin)
            if (all(size(orgin) == size(obj.n)))
                obj.org = orgin;
            else
                if isempty(obj.n)
                    throw(MException('MATLAB:argchk','Define n first'));
                else
                    throw(MException('MATLAB:argchk', ...
                    'Origin must have the same dimensions as the number of points.'));
                end
            end
        end
        
        function obj = setSpacing(obj,spacing)
            if (all(size(spacing) == size(obj.n)))
                obj.spacing = spacing;
            else
                if isempty(obj.n)
                    throw(MException('MATLAB:argchk','Define n first'));
                else
                    throw(MException('MATLAB:argchk', ...
                    'Spacing must have the same dimensions as the number of points.'));
                end
            end
        end
        
        function out = getNumPoints(obj)
            out = prod(obj.n);
        end
        
        function out = okToWrite(obj)
            if (isempty(obj.n) || isempty(obj.org) || isempty(obj.spacing))
                out = false;
            else
                out = true;
            end
        end
        
        function writeFile(obj,filename)
            fid = fopen(filename, 'a'); 
            fprintf(fid,['DATASET STRUCTURED_POINTS\n' ...
                'DIMENSIONS ' num2str(obj.n,'%g ') '\n'...
                'ORIGIN ' num2str(obj.org,'%g ') '\n'...
                'SPACING ' num2str(obj.spacing,'%g ') '\n']);
            fclose(fid);
        end
    end
    
end