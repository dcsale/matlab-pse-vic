function savevtkvector(mesh, X, Y, Z, filename)
%  savevtkvector Save a 3-D vector array in VTK format
%  savevtkvector(X,Y,Z,filename) saves a 3-D vector of any size to
%  filename in VTK format. X, Y and Z should be arrays of the same
%  size, each storing speeds in the a single Cartesian directions.


delX = mesh.h;
delY = mesh.h;
delZ = mesh.h;

if ((size(X) ~= size(Y)) | (size(X) ~= size(Z)))
    fprint('Error: velocity arrays of unequal size \r\n'); return;
end
[nx, ny, nz] = size(mesh.X);
fid = fopen(filename, 'w', 'b');
fprintf(fid, '# vtk DataFile Version 3.0 \r\n');
fprintf(fid, 'Comment goes here \r\n');
fprintf(fid, 'ASCII \r\n');
fprintf(fid, ' \r\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS \r\n');
fprintf(fid, 'DIMENSIONS    %d   %d   %d \r\n', nx, ny, nz);
fprintf(fid, '\r\n');
fprintf(fid, '%s%f%c%f%c%f \r\n', 'ORIGIN ', min(min(min(mesh.X))), ' ', min(min(min(mesh.Y))), ' ', min(min(min(mesh.Z)))); 
fprintf(fid, '%s%f%c%f%c%f \r\n', 'SPACING ', delX, ' ', delY, ' ', delZ); 
fprintf(fid, '\n');
fprintf(fid, 'POINT_DATA   %d \r\n', nx*ny*nz);
fprintf(fid, 'VECTORS vectors float\r\n');
fprintf(fid, '\r\n');
for a=1:nx
    for b=1:ny
        for c=1:nz
            fprintf(fid, '%f ', X(a,b,c));
            fprintf(fid, '%f ', Y(a,b,c));
            fprintf(fid, '%f ', Z(a,b,c));
        end
        fprintf(fid, '\r\n');
    end
end
fclose(fid);

end % function

%%
% function savevtk(array, filename)
% %  savevtk Save a 3-D scalar array in VTK format.
% %  savevtk(array, filename) saves a 3-D array of any size to
% %  filename in VTK format.
%     [nx, ny, nz] = size(array);
%     fid = fopen(filename, 'wt');
%     fprintf(fid, '# vtk DataFile Version 2.0\n');
%     fprintf(fid, 'Comment goes here\n');
%     fprintf(fid, 'ASCII\n');
%     fprintf(fid, '\n');
%     fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
%     fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
%     fprintf(fid, '\n');
%     fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
%     fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
%     fprintf(fid, '\n');
%     fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
%     fprintf(fid, 'SCALARS scalars float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     fprintf(fid, '\n');
%     for a=1:nx
%         for b=1:ny
%             for c=1:nz
%                 fprintf(fid, '%d ', array(a,b,c));
%             end
%             fprintf(fid, '\n');
%         end
%     end
%     fclose(fid);
% return