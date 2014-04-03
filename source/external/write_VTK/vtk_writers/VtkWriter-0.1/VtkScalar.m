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

classdef VtkScalar < handle
    properties (SetAccess=protected)
        comment
        components
        lookuptable
        data
    end
    
    methods
        function obj = VtkScalar(comment,lookuptable)
            obj.comment = comment;
            obj.lookuptable = lookuptable;
        end
        
        function addScalar(obj,varargin)
            start=1;
            if strcmp(varargin{1},'components')
                obj.components = varargin{2};
                start=3;
            end
            
            if (length(varargin)-start~=obj.components && ...
                length(varargin)-start~=0)
                throw(MException('MATLAB:argchk', ...
                    'data dimensions do not match number of components'));
            end 
            if length(varargin{start:end})~=1
                obj.data = [];
                if start == length(varargin)
                    obj.data = reshape(varargin{start},1,numel(varargin{start}));
                else
                    for ii=varargin{start:end}
                        obj.data = [obj.data; reshape(ii,1,numel(ii))];
                    end
                end
            else
                if (size(varargin{start},1)>4)
                    throw(MException('MATLAB:argchk', ...
                    ['Scalars must be 1, 2, 3 or 4 dimensional with '...
                    'each dimension in a row.']));
                end
                obj.components = size(varargin{start},1);
                obj.data = varargin{start};
            end
        end
        
        function writeFile(obj,filename)
            if ~isempty(obj.data)
                fid = fopen(filename, 'a');
                fprintf(fid,['\nSCALARS ' obj.comment ' float ' num2str(obj.components) '\n']);
                fprintf(fid,['LOOKUP_TABLE ' obj.lookuptable '\n']);
                fwrite(fid,obj.data,'float','b');
                fclose(fid);
            end
        end
        
    end
end