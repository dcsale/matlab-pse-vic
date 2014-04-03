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

classdef PointData < handle
    % POINTDATA Visualization Toolkit Point data dataset type.
    % Implements the VTK POINTDATA dataset type. This is a container class
    % for adding dataset attributes to a VTK file.
    %
    % You should not need to call methods of this class! Access is provided
    % in the main VtkWriter class.
    %
    % Currently implemented:
    %   Vectors
    %   Scalars
    %
    % Not Yet Implemented:
    %   Color Scalars
    %   Lookup Table
    %   Normals
    %   Texture_Coordinates
    %   Tensors
    %   Field data
    properties (SetAccess=protected)
        numpoints
        datasets={}
    end
    methods
        function obj = PointData(num)
            obj.numpoints=num;
        end
        
        function addVectorDataset(obj,comment,varargin)
            obj.datasets = [obj.datasets {VtkVector(comment)}];
            obj.datasets{end}.addVector(varargin{:});
        end
        
        function addScalarDataset(obj,comment,lookuptable,varargin)
            obj.datasets = [obj.datasets {VtkScalar(comment,lookuptable)}];
            obj.datasets{end}.addScalar(varargin{:});
        end
        
        function writeHeader(obj,filename)
            fid = fopen(filename, 'a');
            fprintf(fid, ['\nPOINT_DATA ' num2str(obj.numpoints)]);
            fclose(fid);
        end
        
        function writeFile(obj,filename)
            for ii=1:length(obj.datasets)
                obj.datasets{ii}.writeFile(filename);
            end
        end
        
        function clearData(obj)
            obj.datasets={};
        end
    end
end