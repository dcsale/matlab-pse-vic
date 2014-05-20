function plot_particles(xp, wp, MESH, tag)

plot_components = false;

if isempty(wp)
    % can only plot the particle positions
    figure
    set(gcf,'name', tag, 'numbertitle', 'off') 
    
    plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
    axis equal
    axis([MESH.xmin(1) MESH.xmax(1) ...
          MESH.xmin(2) MESH.xmax(2) ...
          MESH.xmin(3) MESH.xmax(3)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    return
end

%% Plot particle weights
wp_mag = sqrt( wp(1,:).^2 + wp(2,:).^2 + wp(3,:).^2 );

figure
set(gcf,'name', tag, 'numbertitle', 'off') 

if plot_components

    subplot(2,2,1)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp(1,:))
    axis equal
    axis([MESH.xmin(1) MESH.xmax(1) ...
          MESH.xmin(2) MESH.xmax(2) ...
          MESH.xmin(3) MESH.xmax(3)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,2)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp(2,:))
    axis equal
    axis([MESH.xmin(1) MESH.xmax(1) ...
          MESH.xmin(2) MESH.xmax(2) ...
          MESH.xmin(3) MESH.xmax(3)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,3)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp(3,:))
    axis equal
    axis([MESH.xmin(1) MESH.xmax(1) ...
          MESH.xmin(2) MESH.xmax(2) ...
          MESH.xmin(3) MESH.xmax(3)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,4)
end

plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp_mag)
axis equal
axis([MESH.xmin(1) MESH.xmax(1) ...
      MESH.xmin(2) MESH.xmax(2) ...
      MESH.xmin(3) MESH.xmax(3)]);
xlabel('x')
ylabel('y')
zlabel('z')

end % plot_particles

