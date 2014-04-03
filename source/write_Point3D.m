function write_Point3D(nPart, xp, wp, up, iter, outputDir, tag)
% Point3D is a simple ASCII text format for 3D particles that have 1 scalar variable. 
% Point3D is not a widely used format but it can be useful for scripts or simple particle 
% simulations that want an easy way to get particle data into VisIt. The file format 
% consists of a comment line that includes the names of the coordinates and the variable 
% name. 
% 
% The variable name is used in VisIt's plot menus. After the comment line, each 
% subsequent line contains the (x,y,z) coordinates and variable value for a point. 
% If you have more than 1 scalar variable to plot on your particles then consider using 
% Silo or Xmdv file formats.  File should have *.3D extension.
% 
% Example file:
% x y z variablename
% 0 0 0 10
% 1 2 3 20
% 2 4 6 30

%% =======================================================================%
% Write particle vorticity
% ========================================================================%
% fid = fopen([outputDir filesep 'particle_vorticity_x_' num2str(iter,'%04.0f') '.3D'],'w','b'); %note: 'b' uses Big-endian ordering
% fprintf(fid, 'xp_x xp_y xp_z wp_x \r\n');
% for n = 1:nPart
%     fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), wp(1,n));
% end
% fclose(fid);
% 
% fid = fopen([outputDir filesep 'particle_vorticity_y_' num2str(iter,'%04.0f') '.3D'],'w','b'); 
% fprintf(fid, 'xp_x xp_y xp_z wp_y \r\n');
% for n = 1:nPart
%     fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), wp(2,n));
% end
% fclose(fid);
% 
% fid = fopen([outputDir filesep 'particle_vorticity_z_' num2str(iter,'%04.0f') '.3D'],'w','b'); 
% fprintf(fid, 'xp_x xp_y xp_z wp_z \r\n');
% for n = 1:nPart
%     fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), wp(3,n));
% end
% fclose(fid);
% 
fid = fopen([outputDir filesep tag '-particle_vorticity_mag_' num2str(iter,'%04.0f') '.3D'],'w','b'); 
fprintf(fid, 'xp_x xp_y xp_z wp_mag \r\n');
for n = 1:nPart
    fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), sqrt(wp(1,n)^2 + wp(2,n)^2 + wp(3,n)^2));
end
fclose(fid);

%% =======================================================================%
% Write particle velocity
% ========================================================================%
% fid = fopen([outputDir filesep 'particle_velocity_x_' num2str(iter,'%04.0f') '.3D'],'w','b'); %note: 'b' uses Big-endian ordering
% fprintf(fid, 'xp_x xp_y xp_z up_x \r\n');
% for n = 1:nPart
%     fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), up(1,n));
% end
% fclose(fid);
% 
% fid = fopen([outputDir filesep 'particle_velocity_y_' num2str(iter,'%04.0f') '.3D'],'w','b'); 
% fprintf(fid, 'xp_x xp_y xp_z up_y \r\n');
% for n = 1:nPart
%     fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), up(2,n));
% end
% fclose(fid);
% 
% fid = fopen([outputDir filesep 'particle_velocity_z_' num2str(iter,'%04.0f') '.3D'],'w','b'); 
% fprintf(fid, 'xp_x xp_y xp_z up_z \r\n');
% for n = 1:nPart
%     fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), up(3,n));
% end
% fclose(fid);

fid = fopen([outputDir filesep tag '-particle_velocity_mag_' num2str(iter,'%04.0f') '.3D'],'w','b'); 
fprintf(fid, 'xp_x xp_y xp_z up_mag \r\n');
for n = 1:nPart
    fprintf(fid, '%8.8f %8.8f %8.8f %8.8f \r\n', xp(1,n), xp(2,n), xp(3,n), sqrt(up(1,n)^2 + up(2,n)^2 + up(3,n)^2));
end
fclose(fid);

end % function write_Point3D
