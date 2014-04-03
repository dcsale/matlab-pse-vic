function plot_scalar(f, MESH, tag)

% convert from NDGRID to MESHGRID format for visualization
xf    = permute(MESH.xf{1}, [2,1,3]);
yf    = permute(MESH.xf{2}, [2,1,3]);
zf    = permute(MESH.xf{3}, [2,1,3]);
f     = permute(         f, [2,1,3]);

% parameters for plotting
nLevels = 11;

%% Field Slices
figure
set(gcf,'name', [tag ' Slices'], 'numbertitle', 'off')          
xyz_slice = [0, 0, 0];

h = slice(xf, yf, zf, f, xyz_slice(1), xyz_slice(2), xyz_slice(3),'cubic');
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
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
drawnow

end % function

