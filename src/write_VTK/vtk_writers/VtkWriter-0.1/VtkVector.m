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

classdef VtkVector < handle
    % VTKVECTOR Visualization Toolkit Vector
    properties (SetAccess=protected)
        comment
        data = []
    end
    
    methods
        function obj = VtkVector(comment)
            % VtkVector Constructor.
            % This is an internal class and it should not be necessary to
            % use the constructor.
            obj.comment = comment;
        end
        
        function addVector(obj,varargin)
            if ~isempty(varargin)
                for ii=1:length(varargin)
                    n = varargin{ii};
                    obj.data = [obj.data; ...
                        reshape(n,1,numel(n))];
                end
            else
                throw(MException('MATLAB:argchk', ...
                    'Cannot add an empty vector'));
            end
        end
        
        function writeFile(obj,filename)
            if ~isempty(obj.data)
                fid = fopen(filename, 'a');
                fprintf(fid,['\nVECTORS ' obj.comment ' float\n']);
                fwrite(fid,obj.data,'float','b');
                fclose(fid);
            end
        end
    end
end