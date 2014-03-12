function write_VTK_v2(varName, mesh, u, v, w, SIM.caseName, SIM.outputDir)

%% =======================================================================%
% Preamble
% ========================================================================%

vtkEndian = 'ieee-le'; % for LittleEndian
% vtkEndian = 'ieee-be'; % for BigEndian

vtkVersion = 'vtk_legacy';
% version = 'vtk_xml';

vtkFileType = 'ImageData'; % for ImageData (.vti) â€” Serial vtkImageData (structured)
% vtkFileType = 'RectilinearGrid'; % for RectilinearGrid (.vtr) â€” Serial vtkRectilinearGrid (structured)
% vtkFileType = 'StructuredGrid'; % for StructuredGrid (.vts) â€” Serial vtkStructuredGrid (structured).

x1 = mesh.x(1);
x2 = mesh.x(end);
y1 = mesh.y(1);
y2 = mesh.y(end);
z1 = mesh.z(1);
z2 = mesh.z(end);
x0 = mesh.x(1);
y0 = mesh.y(1);
z0 = mesh.z(1);
dx = mesh.h;
dy = mesh.h;
dz = mesh.h;

switch vtkVersion
    case 'vtk_legacy'
        fid = fopen([SIM.outputDir filesep SIM.caseName '.vtk'], 'w', vtkEndian);
        
        fprintf(fid, '# vtk DataFile Version 2.0\n');
        fprintf(fid, 'VTK from Matlab\n');
        fprintf(fid, 'ASCII\n');
        
        fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
        fprintf(fid, 'DIMENSIONS %g %g %g\n', numel(mesh.x)-1, numel(mesh.y)-1, numel(mesh.z)-1);
        fprintf(fid, 'ORIGIN %g %g %g\n', x0, y0, z0);
        fprintf(fid, 'SPACING %g %g %g\n', dx, dy, dz);
        
        fprintf(fid, 'POINT_DATA %g\n', numel(mesh.x)*numel(mesh.y)*numel(mesh.z));
        fprintf(fid, 'VECTORS vel_test float\n');
        fprintf(fid, '%g %g %g\n', u, v, w);
        
        fclose(fid);       
         
    case 'vtk_xml'
        switch vtkFileType
            case 'ImageData'
                % Each ImageData piece specifies its extent within the datasetâ€™s whole extent. The points and cells
                % are described implicitly by the extent, origin, and spacing. Note that the origin and spacing are constant across all
                % pieces, so they are specified as attributes of the ImageData XML element as follows.       
                fid = fopen([SIM.outputDir filesep SIM.caseName '.vti'], 'w', vtkEndian);
                fprintf(fid, '<?xml version="1.0"?> \r\n');
                fprintf(fid, '<VTKFile type="%s" version="0.1" byte_order="%s"> \r\n', vtkFileType, vtkEndian);
        %         fprintf(fid, '<VTKFile type="%s" version="0.1"> \r\n', vtkFileType);        
                fprintf(fid, '  <ImageData WholeExtent="%g %g %g %g %g %g" Origin="%g %g %g" Spacing="%g %g %g"> \r\n', 0, numel(mesh.x)-1, 0, numel(mesh.y)-1, 0, numel(mesh.z)-1, x0, y0, z0, dx, dy, dz);
                fprintf(fid, '    <Piece Extent="%g %g %g %g %g %g"> \r\n', 0, numel(mesh.x)-1, 0, numel(mesh.y)-1, 0, numel(mesh.z)-1);

                fprintf(fid,['      <PointData Vectors="' varName '"> \r\n']);
                fprintf(fid,['        <DataArray Name=â€?' varName 'â€? type="Float64" NumberOfComponents="3" format="ascii"> \r\n']);
                fprintf(fid, '        %g %g %g \r\n', u, v, w);
                fprintf(fid, '        </DataArray> \r\n');
                fprintf(fid, '      </PointData> \r\n');

        %         fprintf(fid,['      <PointData Scalars="' varName '"> \r\n']);
        %         fprintf(fid,['        <DataArray Name="' varName '_x" type="Float64" format="ascii"> \r\n']);
        %         fprintf(fid, '        %g \r\n', u);
        %         fprintf(fid, '        </DataArray> \r\n');
        %         fprintf(fid,['        <DataArray Name="' varName '_y" type="Float64" format="ascii"> \r\n']);
        %         fprintf(fid, '        %g \r\n', v);
        %         fprintf(fid, '        </DataArray> \r\n');
        %         fprintf(fid,['        <DataArray Name="' varName '_z" type="Float64" format="ascii"> \r\n']);
        %         fprintf(fid, '        %g \r\n', w);
        %         fprintf(fid, '        </DataArray> \r\n'); 
        %         fprintf(fid, '      </PointData> \r\n');

                fprintf(fid, '      <CellData> \r\n');
                fprintf(fid, '      </CellData> \r\n');

                fprintf(fid, '    </Piece> \r\n');
                fprintf(fid, '  </ImageData> \r\n');
                fprintf(fid, '</VTKFile> \r\n');
                fclose(fid);

            case 'RectilinearGrid'
                % Each RectilinearGrid piece specifies its extent within the datasetï¿½s whole extent. The
                % points are described by the Coordinates element. The cells are described implicitly by the extent.
                fid = fopen([SIM.outputDir filesep SIM.caseName '.vtr'], 'w', vtkEndian);
                fprintf(fid, '<?xml version="1.0"?> \r\n');
                fprintf(fid, '<VTKFile type="%s" version="0.1" byte_order="%s"> \r\n', vtkFileType, vtkEndian);
                fprintf(fid, '  <RectilinearGrid WholeExtent="%g %g %g %g %g %g"> \r\n', x1, x2, y1, y2, z1, z2);
                fprintf(fid, '    <Piece Extent="%g %g %g %g %g %g"> \r\n', x1, x2, y1, y2, z1, z2);

                fprintf(fid,['      <PointData Vectors="' varName '"> \r\n']);
        %         fprintf(fid,['        <DataArray Name=â€?' varName 'â€? type="Float64" NumberOfComponents="3" format="ascii"> \r\n']);
        %         fprintf(fid, '        %g %g %g \r\n', u, v, w);
        %         fprintf(fid, '        </DataArray> \r\n');
                fprintf(fid,['        <DataArray Name=â€?' varName '_xâ€? type="Float64" NumberOfComponents="1" format="ascii"> \r\n']);
                fprintf(fid, '        %g \r\n', u);
                fprintf(fid, '        </DataArray> \r\n');
                fprintf(fid,['        <DataArray Name=â€?' varName '_yâ€? type="Float64" NumberOfComponents="1" format="ascii"> \r\n']);
                fprintf(fid, '        %g \r\n', v);
                fprintf(fid, '        </DataArray> \r\n');
                fprintf(fid,['        <DataArray Name=â€?' varName '_zâ€? type="Float64" NumberOfComponents="1" format="ascii"> \r\n']);
                fprintf(fid, '        %g \r\n', w);
                fprintf(fid, '        </DataArray> \r\n'); 
                fprintf(fid, '      </PointData> \r\n');

                fprintf(fid, '      <CellData> \r\n');
                fprintf(fid, '      </CellData> \r\n');

                % The Coordinates element defines point coordinates for an extent by specifying the ordinate
                % along each axis for each integer value in the extentâ€™s range. It contains three DataArray elements 
                % describing the ordinates along the x-y-z axes, respectively.
                fprintf(fid, '      <Coordinates> \r\n'); 
                fprintf(fid, '        <DataArray Name=â€?mesh_xâ€? type="Float64" NumberOfComponents="1" format="ascii"> \r\n');
                fprintf(fid, '        %g \r\n', mesh.x);
                fprintf(fid, '        </DataArray> \r\n');
                fprintf(fid, '        <DataArray Name=â€?mesh_yâ€? type="Float64" NumberOfComponents="1" format="ascii"> \r\n');
                fprintf(fid, '        %g \r\n', mesh.y);
                fprintf(fid, '        </DataArray> \r\n');
                fprintf(fid, '        <DataArray Name=â€?mesh_zâ€? type="Float64" NumberOfComponents="1" format="ascii"> \r\n');
                fprintf(fid, '        %g \r\n', mesh.z);
                fprintf(fid, '        </DataArray> \r\n');
                fprintf(fid, '      </Coordinates> \r\n');

                fprintf(fid, '    </Piece> \r\n');
                fprintf(fid, '  </RectilinearGrid> \r\n');
                fprintf(fid, '</VTKFile> \r\n');
                fclose(fid);

            case 'StructuredGrid'
                % Each StructuredGrid piece specifies its extent within the datasetâ€™s whole extent. The points
                % are described explicitly by the Points element. The cells are described implicitly by the extent.
                fid = fopen([SIM.outputDir filesep SIM.caseName '.vts'], 'w', vtkEndian);
                fprintf(fid, '<?xml version="1.0"?> \r\n');
                fprintf(fid,['<VTKFile type="' vtkFileType '" version="0.1" byte_order="' vtkEndian '"> \r\n']);
                fprintf(fid, '  <StructuredGrid WholeExtent="%g %g %g %g %g %g"> \r\n', x1, x2, y1, y2, z1, z2);
                fprintf(fid, '    <Piece Extent="%g %g %g %g %g %g"> \r\n', x1, x2, y1, y2, z1, z2);
                fprintf(fid,['      <PointData Scalars="' Temperature '"> \r\n']);        
                fprintf(fid, '      </PointData> \r\n');
                fprintf(fid, '      <CellData> \r\n');
                fprintf(fid, '      </CellData> \r\n');
                fprintf(fid, '      <Points> \r\n');

                fprintf(fid, '      </Points> \r\n');
                fprintf(fid, '    </Piece> \r\n');
                fprintf(fid, '  </StructuredGrid> \r\n');
                fprintf(fid, '</VTKFile> \r\n');
                fclose(fid);

            otherwise
                error('[ERROR] Could not recognize input for variable: vtkFileType')

        end % switch vtkFileType
end

end % function
