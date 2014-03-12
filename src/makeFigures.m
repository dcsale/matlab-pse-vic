function makeFigures(mesh, ring, xp, wp, up, uf_x, uf_y, uf_z, wf_x, wf_y, wf_z, DEBUG_LVL)

uf_mag = sqrt( uf_x.^2 + uf_y.^2 + uf_z.^2 );
wf_mag = sqrt( wf_x.^2 + wf_y.^2 + wf_z.^2 );
wp_mag = sqrt( wp(1,:).^2 + wp(2,:).^2 + wp(3,:).^2 );
up_mag = sqrt( up(1,:).^2 + up(2,:).^2 + up(3,:).^2 );

% Visualization

%% Vorticity Field
if DEBUG_LVL > 1
    figure
    set(gcf,'name', 'Vorticity Field','numbertitle','off')
    subplot(2,2,1)
    h = slice(mesh.X, mesh.Y, mesh.Z, wf_x, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,2)
    h = slice(mesh.X, mesh.Y, mesh.Z, wf_y, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,3)
    h = slice(mesh.X, mesh.Y, mesh.Z, wf_z, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,4)
    h = slice(mesh.X, mesh.Y, mesh.Z, wf_mag, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow

    %% Vorticity Iso-Surfaces
    figure
    set(gcf,'name', 'Vorticity Isosurfaces','numbertitle','off')
    % isosurface(X,Y,Z,W,0)
    isosurface(mesh.X,mesh.Y,mesh.Z,wf_mag, 0.25 * max(max(max(wf_mag))))
    isosurface(mesh.X,mesh.Y,mesh.Z,wf_mag, 0.50 * max(max(max(wf_mag))))
    isosurface(mesh.X,mesh.Y,mesh.Z,wf_mag, 0.75 * max(max(max(wf_mag))))
    alpha(0.25)
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')

end

if DEBUG_LVL > 2    
    %% Velocity Field
    figure
    set(gcf,'name', 'Velocity Field','numbertitle','off')
    subplot(2,2,1)
    h = slice(mesh.X, mesh.Y, mesh.Z, uf_x, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    % plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,2)
    h = slice(mesh.X, mesh.Y, mesh.Z, uf_y, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    % plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,3)
    h = slice(mesh.X, mesh.Y, mesh.Z, uf_z, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    % plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,4)
    h = slice(mesh.X, mesh.Y, mesh.Z, uf_mag, [ring.center_x], [ring.center_y], [ring.center_z],'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    % plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')

    %% Particle Vorticity
    figure
    set(gcf,'name', 'Particle Vorticity','numbertitle','off')
    subplot(2,2,1)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp(1,:))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,2)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp(2,:))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,3)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp(3,:))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,4)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp_mag)
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')

    % Particle Velocity
    figure(5)
    clf
    set(gcf,'name', 'Particle Velocity','numbertitle','off')
    subplot(2,2,1)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(1,:))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow

    subplot(2,2,2)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(2,:))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow

    subplot(2,2,3)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(3,:))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow

    subplot(2,2,4)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up_mag)
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow
end