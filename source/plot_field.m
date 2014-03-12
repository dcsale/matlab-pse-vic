function plot_field(field, MESH, tag)


%% some pre-processing of the field data
f_x   = field{1};
f_y   = field{2};
f_z   = field{3};
f_mag = sqrt( f_x.^2 + f_y.^2 + f_z.^2 );

% convert from NDGRID to MESHGRID format for visualization
MESH.xf = permute(MESH.xf, [2,1,3]);
MESH.yf = permute(MESH.yf, [2,1,3]);
MESH.zf = permute(MESH.zf, [2,1,3]);
f_x     = permute(    f_x, [2,1,3]);
f_y     = permute(    f_y, [2,1,3]);
f_z     = permute(    f_z, [2,1,3]);
f_mag   = permute(  f_mag, [2,1,3]);

% parameters for plotting
nLevels = 11;
%         cmap = colormap(bipolar(2*nLevels, 0.51));

%% Field Slices
figure
set(gcf,'name', [tag ' Slices'], 'numbertitle', 'off')      
[f_mag_max, imax] = max(f_mag(:));        
%         xyz_slice = [0, 0, 0];
xyz_slice = [MESH.xf(imax), MESH.yf(imax), MESH.zf(imax)];
subplot(2,2,1)
h = slice(MESH.xf, MESH.yf, MESH.zf, f_x, xyz_slice(1), xyz_slice(2), xyz_slice(3),'cubic');
set(h,'FaceColor',          'interp', ...
      'EdgeColor',          'none', ...
      'DiffuseStrength',    0.8)
cb = colorbar;
colormap(bipolar(2*nLevels, 0.51));
cabs = max(abs(get(gca, 'CLim')));
set(gca, 'CLim',  [-cabs cabs]);
set(cb,  'YLim',  [-cabs cabs]);
set(cb,  'YTick', consolidator([linspace(-cabs,    0, nLevels/2 + 1), ...
                                linspace(    0, cabs, nLevels/2 + 1)]));
hold on
% plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,2,2)
h = slice(MESH.xf, MESH.yf, MESH.zf, f_y, xyz_slice(1), xyz_slice(2), xyz_slice(3),'cubic'); 
set(h,'FaceColor',          'interp', ...
      'EdgeColor',          'none', ...
      'DiffuseStrength',    0.8)
cb = colorbar;
colormap(bipolar(2*nLevels, 0.51));
cabs = max(abs(get(gca, 'CLim')));
set(gca, 'CLim',  [-cabs cabs]);
set(cb,  'YLim',  [-cabs cabs]);
set(cb,  'YTick', consolidator([linspace(-cabs,    0, nLevels/2 + 1), ...
                                linspace(    0, cabs, nLevels/2 + 1)]));
hold on
% plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
subplot(2,2,3)
h = slice(MESH.xf, MESH.yf, MESH.zf, f_z, xyz_slice(1), xyz_slice(2), xyz_slice(3),'cubic');
set(h,'FaceColor',          'interp', ...
      'EdgeColor',          'none', ...
      'DiffuseStrength',    0.8)
cb = colorbar;
colormap(bipolar(2*nLevels, 0.51));
cabs = max(abs(get(gca, 'CLim')));
set(gca, 'CLim',  [-cabs cabs]);
set(cb,  'YLim',  [-cabs cabs]);
set(cb,  'YTick', consolidator([linspace(-cabs,    0, nLevels/2 + 1), ...
                                linspace(    0, cabs, nLevels/2 + 1)]));
hold on
% plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
subplot(2,2,4)
h = slice(MESH.xf, MESH.yf, MESH.zf, f_mag, xyz_slice(1), xyz_slice(2), xyz_slice(3),'cubic');
set(h,'FaceColor',          'interp', ...
      'EdgeColor',          'none', ...
      'DiffuseStrength',    0.8)
cb = colorbar;
colormap(bipolar(2*nLevels, 0.51));
cabs = max(abs(get(gca, 'CLim')));
set(gca, 'CLim',  [-cabs cabs]);
set(cb,  'YTick', linspace(0, cabs, nLevels+1));
set(cb,  'YLim',  [0 cabs]);
hold on
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
drawnow    

%% Iso-Surfaces
figure
set(gcf,'name', [tag ' Isosurfaces'],'numbertitle','off')
clf
isosurface(MESH.xf, MESH.yf, MESH.zf, f_mag, 0.25 * f_mag_max)
isosurface(MESH.xf, MESH.yf, MESH.zf, f_mag, 0.50 * f_mag_max)
isosurface(MESH.xf, MESH.yf, MESH.zf, f_mag, 0.75 * f_mag_max)
hold on
%         plot3(xp(1,:), ...
%               xp(2,:), ...
%               xp(3,:), ...
%               '.g');
alpha(0.25)
axis equal
axis([MESH.xmin(1) MESH.xmax(1) ...
      MESH.xmin(2) MESH.xmax(2) ...
      MESH.xmin(3) MESH.xmax(3)]);
xlabel('x')
ylabel('y')
zlabel('z')
   
end

