function write_vector_field_VTK(varName, mesh, u, v, w, time, format, prec, vtkFile)

switch prec
    case 'single'
        Points  = [single(mesh.X(:)) single(mesh.Y(:)) single(mesh.Z(:))];     % grid coordinates
        Vectors = [single(     u(:)) single(     v(:)) single(     w(:))];     % vector at every grid point
    case 'double'
        Points  = [double(mesh.X(:)) double(mesh.Y(:)) double(mesh.Z(:))];     % grid coordinates
        Vectors = [double(     u(:)) double(     v(:)) double(     w(:))];     % vector at every grid point
end

% open file
fid = fopen(vtkFile,'w','ieee-le');

%% == Header ==
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'ASCII\n');


%% mesh data

% == Structured Grid == 
% fprintf(fid, 'DATASET STRUCTURED_GRID\n');
% fprintf(fid, ['DIMENSIONS ' num2str(size(mesh.X,1)) ' ' num2str(size(mesh.X,2)) ' ' num2str(size(mesh.X,3)) '\n']);
% fprintf(fid, ['POINTS ' num2str(numel(mesh.X)) ' double\n']);
% for i = 1:size(Points,1)
%     fprintf(fid, '%f %f %f \n', Points(i,1), Points(i,2), Points(i,3));
% end

% == Structured Points ==
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
fprintf(fid, ['DIMENSIONS ' num2str(size(mesh.X,1)) ' ' num2str(size(mesh.X,2)) ' ' num2str(size(mesh.X,3)) '\n']);
fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', min(min(min(mesh.X))), ' ', min(min(min(mesh.Y))), ' ', min(min(min(mesh.Z)))); 
fprintf(fid, '%s%f%c%f%c%f\n', 'SPACING ', mesh.h, ' ', mesh.h, ' ', mesh.h); 

%% field data
% fprintf(fid, ['FIELD ' varName ' 3\n']);
% fprintf(fid, 'TIME 1 1 double\n');
% fprintf(fid, '666\n');
% fprintf(fid, 'CYCLE 1 1 int\n');
% fprintf(fid, '9000\n');
% fprintf(fid, ['vel_test 3 ' num2str(numel(mesh.X)) ' double\n']);
fprintf(fid, 'VECTORS vel_test double\n\n');
for i = 1:size(Vectors,1)
    fprintf(fid, '%f %f %f \n', Vectors(i,1), Vectors(i,2), Vectors(i,3));
end

%%
fclose(fid);



% 
% % fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', min(min(min(mesh.X))), ' ', min(min(min(mesh.Y))), ' ', min(min(min(mesh.Z)))); 
% % fprintf(fid, '%s%f%c%f%c%f\n', 'SPACING ', mesh.h, ' ', mesh.h, ' ', mesh.h); 
% 
% % fprintf(fid, 'DATASET STRUCTURED_GRID\n');
% 
% 
% fprintf(fid, '%s%d\n', 'POINT_DATA ', numel(mesh.X));
% fclose(fid);
% 
% %append binary x,y,z data
% fid = fopen(vtkfile, 'a');
% %append binary x,y,z data
% fwrite(fid, [reshape(mesh.X,1,nr_of_elements);  reshape(mesh.Y,1,nr_of_elements); reshape(mesh.Z,1,nr_of_elements)],'float','b');
% 
% %append another ASCII sub header
% fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
% fprintf(fid, 'VECTORS vectors double\n');
% 
% %append binary u,v,w data
% fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b');
% 
% %append some scalar data
% fprintf(fid, '\nSCALARS vector_mag double\n'); %ASCII header
% fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
% fprintf(fid, reshape(sqrt(u.^2+v.^2+w.^2),1,nr_of_elements),'double','b');

% fclose(fid);

% 
% 
% 
% 
% 
% 
% %append another ASCII sub header
% fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
% fprintf(fid, 'VECTORS field float\n');
% 
% %append binary u,v,w data
% fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b');
% 
% %append another ASCII sub header
% fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
% fprintf(fid, 'VECTORS velocity_vectors float\n');
% 
% %append some scalar data
% % fprintf(fid, '\nSCALARS vector_mag float\n'); %ASCII header
% % fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
% % fwrite (fid, reshape(sqrt(u.^2+v.^2+w.^2),1,nr_of_elements),'float','b'); %binary data
%  
% % append some scalar data
% varname = 'vector_mag';
% varInfo = whos('u');
% tp = varInfo.class;
% if( strcmp(tp, 'uint8') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' unsigned_char']);    
% elseif( strcmp(tp, 'int8') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' char']);    
% elseif( strcmp(tp, 'uint16') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' unsigned_short']);    
% elseif( strcmp(tp, 'int16') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' short']);    
% elseif( strcmp(tp, 'uint32') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' unsigned_int']);    
% elseif( strcmp(tp, 'int32') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' int']);    
% elseif( strcmp(tp, 'single') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' float']);    
% elseif( strcmp(tp, 'double') > 0 )
%     fprintf(fid, '%s\n', ['SCALARS ' varname ' double']);   
% end
% fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
% fwrite (fid, reshape(sqrt(u.^2 + v.^2 + w.^2),1,nr_of_elements),'float','b'); %binary data
% 
% fclose(fid);