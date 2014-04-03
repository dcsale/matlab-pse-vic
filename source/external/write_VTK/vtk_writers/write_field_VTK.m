function write_field_VTK(varName, mesh, u, v, w, time, format, prec, vtkFile)
% An example file that shows how to write a STRUCTURED VTK-XML file from
% matlab, both in ASCII and Binary format, with scalar and vector data.
%
% Boris Kaus, ETH Zurich, Feb 26 2008.

% Transform the arrays in vector-format, and change them to the appropriate
% precision. This example file assumes that data comes in single precision
%
% Double precision requires changing 
switch prec
    case 'single'
        Points  = [single(mesh.X(:)) single(mesh.Y(:)) single(mesh.Z(:))];     % grid coordinates
        Vectors = [single(     u(:)) single(     v(:)) single(     w(:))];     % vector at every grid point
    case 'double'
        Points  = [double(mesh.X(:)) double(mesh.Y(:)) double(mesh.Z(:))];     % grid coordinates
        Vectors = [double(     u(:)) double(     v(:)) double(     w(:))];     % vector at every grid point
end

% Definitions and initialization
sizeof_Float32 = 4;      
sizeof_Float64 = 4;     
sizeof_UInt32  = 4; 
sizeof_UInt64  = 4; 
offset         = 0;      % Initial offset

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fid = fopen(vtkFile, 'w', 'ieee-le');	% 'ieee-be' = BigEndian, 'ieee-le' = LittleEndian
fprintf(fid,'<?xml version="1.0"?> \n');
% fprintf(fid,'<DataSet timestep"%f">\n', time);
fprintf(fid,'<VTKFile type="StructuredGrid" version=”0.1” byte_order=”LittleEndian”>\n');
fprintf(fid,'  <StructuredGrid  WholeExtent="%i %i %i %i %i %i">\n', [0 size(mesh.X,1)-1 0 size(mesh.X,2)-1 0 size(mesh.X,3)-1]);
fprintf(fid,'  <Piece Extent="%i %i %i %i %i %i">\n',                [0 size(mesh.X,1)-1 0 size(mesh.X,2)-1 0 size(mesh.X,3)-1]);

fprintf(fid,['    <PointData Vectors="' varName '">\n']);

%--------------------------------------------------------------------------
% write vector field
switch format
    case 'ASCII'
        fprintf(fid,['      <DataArray Name="' varName '" NumberOfComponents="3" format="ascii" type="Float64">\n']);
        for i = 1:size(Vectors,1)
            fprintf(fid, '%g %g %g \n', Vectors(i,:));
        end
    case 'BINARY'
        fprintf(fid,['      <DataArray Name="' varName '" NumberOfComponents="3" format="appended" type="Float64" offset="%i">\n'],int64(offset));
        offset = offset + length(Vectors(:))*sizeof_Float64 + 1*sizeof_UInt64;  % update the offset       
end
fprintf(fid,'      </DataArray>\n');
% -----------------------
fprintf(fid,'    </PointData>\n');

fprintf(fid,'    <Celldata>\n');
fprintf(fid,'    </Celldata>\n');

%--------------------------------------------------------------------------
% write mesh coordinates
fprintf(fid,'    <Points>\n');
switch format
    case 'ASCII'
        fprintf(fid,'      <DataArray Name="Mesh" NumberOfComponents="3" format="ascii" type="Float64">\n');
        for i = 1:size(Points,1)
            fprintf(fid, '%g %g %g \n', Points(i,:));
        end
    case 'BINARY'
        fprintf(fid,'      <DataArray Name="Mesh" NumberOfComponents="3" format="appended" offset="%i" type="Float64">\n',int64(offset));
        offset = offset + length(Points(:))*sizeof_Float64 + 1*sizeof_UInt64;   % update the offset      
end
fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');
%--------------------------------------------------------------------------

fprintf(fid,'  </Piece> \n');
fprintf(fid,'  </StructuredGrid> \n');

if strcmp(format, 'BINARY')
    % Append binary data in raw format: the order in which data arrays are
    % added should be the same as how they are defined above
    fprintf(fid,'  <AppendedData encoding="raw"> \n');
    fprintf(fid,'_');

    % Add Vectors in binary format
    fwrite(fid, uint32(length(Vectors(:))*sizeof_Float64), 'uint64');
    fwrite(fid, Vectors.', 'float64');

    % Add Points in binary format
    fwrite(fid, uint32(length(Points(:))*sizeof_Float64), 'uint64');
    fwrite(fid, Points.', 'float64');

    fprintf(fid,'  </AppendedData> \n');
end

fprintf(fid,'</VTKFile>\n');
fprintf(fid,'</DataSet>\n');
fclose(fid);
