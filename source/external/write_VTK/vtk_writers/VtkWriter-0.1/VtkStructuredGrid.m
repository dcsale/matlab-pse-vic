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

classdef VtkStructuredGrid < handle
    properties (SetAccess=private)
        n
        points
    end
    methods
        function obj = VtkStructuredGrid(varargin)
            if ~isempty(varargin)
                obj.setWithMeshgrid(varargin{:})
            end
        end
        
        function obj = setGridSize(obj,n)
            
            if (all(size(n)==[1 1]) ||all(size(n)==[1 2]) || all(size(n)==[1 3]))
                obj.n = n;
            else
                throw(MException('MATLAB:argchk',...
                    'Structured Grid only supports 1D, 2D and 3D.'));
            end
        end
        
        function obj = setWithMeshgrid(obj,varargin)
            switch length(varargin)
                case 2
                    [X Y] = meshgrid(varargin{:});
                    obj.points = [reshape(X,1,numel(X)); ...
                        reshape(Y,1,numel(Y))];
                case 3
                    [X Y Z] = meshgrid(varargin{:});
                    obj.points = [reshape(X,1,numel(X)); ...
                        reshape(Y,1,numel(Y)); ...
                        reshape(Z,1,numel(Z))];
                otherwise
                    throw(MException('MATLAB:argchk',...
                        'Meshgrid generation of grid is only supported for 2D and 3D.'));
            end
            obj.n = size(X);
            
        end
        
        function obj = setGrid(obj,varargin)
            switch length(varargin)
                case 1
                    if numel(varargin{1}) ~= prod(obj.n)
                        throw(MException('MATLAB:argchk', ...
                            'Incorrect dimensions of x matrix'));
                    else
                        obj.points = reshape(varargin{1},1,numel(varargin{1}));
                    end
                case 2
                    if numel(varargin{1}) ~= prod(obj.n)
                        throw(MException('MATLAB:argchk', ...
                            'Incorrect dimensions of x matrix'));
                    elseif numel(varargin{2}) ~= prod(obj.n)
                        throw(MException('MATLAB:argchk', ...
                            'Incorrect dimensions of y matrix'));
                    else
                        obj.points = [reshape(varargin{1},1,numel(varargin{1})); ...
                            reshape(varargin{2},1,numel(varargin{2}))];
                    end
                case 3
                    if numel(varargin{1}) ~= prod(obj.n)
                        throw(MException('MATLAB:argchk', ...
                            'Incorrect dimensions of x matrix'));
                    elseif numel(varargin{2}) ~= prod(obj.n)
                        throw(MException('MATLAB:argchk', ...
                            'Incorrect dimensions of y matrix'));
                    elseif numel(varargin{3}) ~= prod(obj.n)
                        throw(MException('MATLAB:argchk', ...
                            'Incorrect dimensions of z matrix'));
                    else
                        obj.points = [reshape(varargin{1},1,numel(varargin{1})); ...
                            reshape(varargin{2},1,numel(varargin{2})); ...
                            reshape(varargin{3},1,numel(varargin{3}))];
                    end
                otherwise
                    throw(MException('MATLAB:argchk',...
                        'Structured Grid is only implemented for 2D and 3D.'));
            end
        end
        
        function out = getNumPoints(obj)
            out = prod(obj.n);
        end
        
        function out = okToWrite(obj)
            if (isempty(obj.n) || isempty(obj.points))
                out = false;
            else
                out = true;
            end
        end
        
        function writeFile(obj,filename)
            if (length(obj.n)==2)
                str = [num2str(obj.n(1)) ' ' num2str(obj.n(2)) ' 1'];
            else
                str = num2str(obj.n,'%i ');
            end
            fid = fopen(filename, 'a');
            fprintf(fid,['DATASET STRUCTURED_GRID\n' ...
                'DIMENSIONS '  str ' \n'...
                'POINTS ' num2str(prod(obj.n),'%i') ' float\n']);
            fwrite(fid,obj.points,'float','b');
            fclose(fid);
        end
        
        function clearGrid(obj)
            obj.points=[];
        end
    end
    
end