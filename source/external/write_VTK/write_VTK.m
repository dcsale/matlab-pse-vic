function [] = write_VTK(X,Y,Z,vol,vtkfile,varname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: writeVTK(vol,vtkfile)
%
%   vol:     The 3D matrix to be saved to file
%   vtkfile: The output filename (string)
%   notes:   Only writes binary STRUCTURED_POINTS
%  
% Erik Vidholm 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions
volinfo = whos('vol');
sz      = volinfo.size;
dimX    = sz(1); 
dimY    = sz(2); 
if length(sz) == 3 
    dimZ = sz(3);
else
    dimZ = 1;
end

% spacing
delX = X(1,2,1) - X(1,1,1);
delY = Y(2,1,1) - Y(1,1,1);
delZ = Z(1,1,2) - Z(1,1,1);

% open file (OBS! big endian format)
fid = fopen(vtkfile,'w','b');

% write header
fprintf(fid, '%s\n', '# vtk DataFile Version 3.0');
fprintf(fid, '%s\n', 'created by write_VTK (Matlab implementation by Erik Vidholm 2006, modified by Danny Sale 2013)');
fprintf(fid, '%s\n', 'BINARY');  
fprintf(fid, '%s\n', 'DATASET STRUCTURED_POINTS');  
fprintf(fid, '%s%d%c%d%c%d\n', 'DIMENSIONS ', dimX, ' ', dimY, ' ', dimZ);
fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', min(min(min(X))), ' ', min(min(min(Y))), ' ', min(min(min(Z)))); 
fprintf(fid, '%s%f%c%f%c%f\n', 'SPACING ', delX, ' ', delY, ' ', delZ); 
fprintf(fid, '%s%d\n', 'POINT_DATA ', dimX*dimY*dimZ);

tp = volinfo.class;
if( strcmp(tp, 'uint8') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' unsigned_char']);
    
elseif( strcmp(tp, 'int8') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' char']);
    
elseif( strcmp(tp, 'uint16') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' unsigned_short']);
    
elseif( strcmp(tp, 'int16') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' short']);
    
elseif( strcmp(tp, 'uint32') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' unsigned_int']);
    
elseif( strcmp(tp, 'int32') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' int']);
    
elseif( strcmp(tp, 'single') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' float']);
    
elseif( strcmp(tp, 'double') > 0 )
    fprintf(fid, '%s\n', ['SCALARS ' varname ' double']);
    
end

fprintf(fid, '%s\n', 'LOOKUP_TABLE default');

% write data as binary
fwrite(fid,vol,tp);

% close file
fclose(fid);

% %%
% %Output file name
% filename=('TestFile.vtk');
% 
% %load the MATLAB 3D flow example
% load wind
% tic
% 
% nr_of_elements=numel(x);
% fid = fopen(filename, 'w'); 
% 
% %ASCII file header
% fprintf(fid, '# vtk DataFile Version 3.0\n');
% fprintf(fid, 'VTK from Matlab\n');
% fprintf(fid, 'BINARY\n\n');
% fprintf(fid, 'DATASET STRUCTURED_GRID\n');
% fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
% fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
% fclose(fid);
% 
% %append binary x,y,z data
% fid = fopen(filename, 'a'); 
% fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');
% 
% %append another ASCII sub header
% fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
% fprintf(fid, 'VECTORS velocity_vectors float\n');
% 
% %append binary u,v,w data
% fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b');
% 
% %append another binary u,v,w data set
% fprintf(fid, '\nVECTORS another_vector_set float\n'); %ASCII header
% fwrite(fid, [reshape(u*10,1,nr_of_elements);  reshape(v*2,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b'); %binary data
% 
% %append some scalar data
% fprintf(fid, '\nSCALARS EinLustigerSkalar float\n'); %ASCII header
% fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
% fwrite (fid, reshape(sqrt(u.^2+v.^2+w.^2),1,nr_of_elements),'float','b'); %binary data
% 
% fclose(fid);
% toc