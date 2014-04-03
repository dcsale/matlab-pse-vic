% VTKWRITER Visualization Toolkit .VTK writer.
% Implements Vtk DataFile Version 3.0. Stores data as floats in binary.
% See http://www.vtk.org/VTK/img/file-formats.pdf for the definition of the standard.
%
% VTKWRITER(filename, comment, type, ...)
% filename - Name of the resulting .vtk file.
% comment - .vtk file title as a string. Must be less than 256 characters.
% type - VtkDataset for this writer.
%
% Additional arguments depend on the type of dataset being used.
%
% Example:
% Creates a .vtk file for the Matlab wind sample.
%   load wind;
%   writer = VtkWriter('wind.vtk','wind example from matlab',VtkDataset.Structured_Grid);
%   writer.grid.setGridSize(size(x));
%   writer.grid.setGrid(x,y,z);
%   writer.addPointData();
%   writer.addVectorToDataset('u_v_w',u,v,w);
%   writer.addScalarToDataset('mag','default','components',1,sqrt(u.^2+v.^2+w.^2))
%   writer.write()
%
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

classdef VtkWriter < handle
    % VTKWRITER Visualization Toolkit .VTK writer.
    % Implements Vtk DataFile Version 3.0. Stores data as floats in binary.
    % See http://www.vtk.org/VTK/img/file-formats.pdf for the definition of
    % the standard.
    %
    % VTKWRITER(filename, comment, type, ...)
    % filename - Name of the resulting .vtk file.
    % comment - .vtk file title as a string. Must be less than 256 characters.
    % type - VtkDataset for this writer.
    %
    % Additional arguments depend on the type of dataset being used.
    %
    % Example:
    % Creates a .vtk file for the Matlab wind sample.
    %   load wind;
    %   writer = VtkWriter('wind.vtk','wind example from matlab',VtkDataset.Structured_Grid);
    %   writer.grid.setGridSize(size(x));
    %   writer.grid.setGrid(x,y,z);
    %   writer.addPointData();
    %   writer.addVectorToDataset('u_v_w',u,v,w);
    %   writer.addScalarToDataset('mag','default','components',1,sqrt(u.^2+v.^2+w.^2))
    %   writer.write()
    
    properties (SetAccess=protected)
        filename='default.vtk';  % Name of the output .vtk file.
        comment='VTK file written by VtkWriter for Matlab'; % .vtk file title as a string. Must be less than 256 characters.
        grid;  % geometry/topology of the .vtk file.
        data; % dataset object (either a PointData or CellData class).
    end
    
    properties (SetAccess=private)
        type % Vtk Dataset type see VtkDataset for valid inputs.
    end
    
    methods
        function obj = VtkWriter(filename,comment,type,varargin)
            % VTKWRITER(filename, comment, type, ...)
            % filename - Name of the resulting .vtk file.
            % comment - .vtk file title as a string. Must be less than 256 characters.
            % type - VtkDataset for this writer.
            %
            % Variable arguments:
            %   VtkDataset.Structured_Grid:
            %       Varable arguments are passed to MESHGRID() for grid
            %       generation.
            %
            %   VtkDataset.Structured_Points:
            %       VTKWRITER(filename,comment,VtkDataset.Structured_Points, n, origin, spacing)
            %       n - number of points. eg. [10 10 21]
            %       origin - start points. eg. [0 0 -1]
            %       spacing - spacing between points. eg. [1 1 0.1]
            if isempty(filename)
                throw(MException('MATLAB:argchk','Invalid Filename'));
            end
            obj.filename=filename;
            if (length(comment) > 255)
                throw(MException('MATLAB:argchk','The comment must be less than 255 characters in length'));
            else
                obj.comment=comment;
            end
            if isa(type,'VtkDataset')
                obj.type=type;
            else
                throw(MException('MATLAB:argchk','The type was not a member of the VtkDataset enumeration'));
            end
            switch type
                case VtkDataset.Structured_Points
                    obj.grid = VtkStructuredPoints(varargin{:});
                case VtkDataset.Structured_Grid
                    obj.grid = VtkStructuredGrid(varargin{:});
                otherwise
                    throw(MException('MATLAB:argchk',[char(type) ': Not Yet Implemented.']));
            end
        end
        
        function setDim(obj,dims)
            % SETDIM(dims)
            if ismember(obj.type,[VtkDataset.Polydata VtkDataset.Unstructured_Grid VtkDataset.Field])
                throw(MException('Invalid:fncall',['Do not use setDim when using a ' char(obj.type) ' VTK Dataset type.']))
            end
            obj.grid.setDims(dims);
        end
        
        function addPointData(obj)
            obj.data = PointData(obj.grid.getNumPoints());
        end
        
        function addVectorToDataset(obj,comment,varargin)
            % addVectorToDataset(comment, ...)
            % Add a vector attribute to the VTK file
            %
            % comment - data name for the vector. No spaces in this comment.
            % 
            % Varable arguments: each component of the vector is given as a separate argument.
            %
            % Example:
            %   writer.addVectorToDataset('att_name',u,v,w);
            obj.data.addVectorDataset(comment,varargin{:});
        end
        
        function addScalarToDataset(obj,comment,lookuptable,varargin)
            % addScalarToDataset(comment,lookuptable,...)
            % Add a scalar attribute to the VTK file.
            %
            % comment - data name for the scalar. No spaces in this comment.
            % lookuptable - Specifies a lookup table to use, typically 'default'.
            %
            % Varable arguments:
            %   addScalarToDataset(comment,lookuptable,numComp,component1,...)
            %   numComp - number of components in the scalar attribute.
            %   further arguments should contain each of the scalar componenets.
            obj.data.addScalarDataset(comment,lookuptable,varargin{:});
        end
        
        function write(obj)
            % write()
            % Write the VTK file at once as currently defined.
            obj.writeHeader()
            obj.data.writeFile(obj.filename);
        end
        
        function writeHeader(obj)
            % writeHeader()
            % Write the header for the VTK file. Useful if dealing with
            % large datasets that will not fit inside the memory space at
            % once.
            fid = fopen(obj.filename, 'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, [obj.comment '\n']);
            fprintf(fid, 'BINARY\n');
            fclose(fid);
            obj.grid.writeFile(obj.filename);
            obj.writeDatasetHeader();
        end
        
        function writeDatasetHeader(obj)
            % Internal function. Do not use. Use writeHeader().
            obj.data.writeHeader(obj.filename);
        end
        
        function writeAndClearData(obj)
            % writeAndClearData()
            % Write the current data to the VTK file and then remove it
            % from the VTK object. You must make a call to writeHeader()
            % before using making calls to this function (otherwise the VTK
            % file will be useless).
            %
            % Usage:
            %   **Initialise writer**
            %
            %   ** Write VTK file header once **
            %   writer.writeHeader()
            %
            %   for each data attribute
            %       ** Load data **
            %
            %       ** Add data to dataset **
            %       eg: writer.addScalarToDataset(...)
            %
            %       ** write data attribute to VTK file **
            %       writer.writeAndClearData()
            %   end
            
            obj.data.writeFile(obj.filename);
            obj.data.clearData();
        end
    end
end