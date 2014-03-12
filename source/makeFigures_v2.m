function makeFigures_v2(MESH, LL, xp, ap, up, uf_x, uf_y, uf_z, wf_x, wf_y, wf_z, DEBUG_LVL)

uf_mag = sqrt( uf_x.^2 + uf_y.^2 + uf_z.^2 );
wf_mag = sqrt( wf_x.^2 + wf_y.^2 + wf_z.^2 );
up_mag = sqrt( up(1,:).^2 + up(2,:).^2 + up(3,:).^2 );
ap_mag = sqrt( ap(1,:).^2 + ap(2,:).^2 + ap(3,:).^2 );


% Visualization

%% Vorticity Field
if DEBUG_LVL > 1
    figure
    set(gcf,'name', 'Vorticity Field','numbertitle','off')
    subplot(2,2,1)
        h = slice(MESH.X, MESH.Y, MESH.Z, wf_x, 0, 0, LL.HUB_HT, 'cubic');
        set(h,'FaceColor',          'interp', ...
              'EdgeColor',          'none', ...
              'DiffuseStrength',    0.8)
        colorbar
        hold on
        plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
        axis equal
        axis([min(MESH.x) max(MESH.x) ...
              min(MESH.y) max(MESH.y) ...
              min(MESH.z) max(MESH.z)]);
        xlabel('x')
        ylabel('y')
        zlabel('z')
    subplot(2,2,2)
        h = slice(MESH.X, MESH.Y, MESH.Z, wf_y, 0, 0, LL.HUB_HT, 'cubic');
        set(h,'FaceColor',          'interp', ...
              'EdgeColor',          'none', ...
              'DiffuseStrength',    0.8)
        colorbar
        hold on
        plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
        axis equal
        axis([min(MESH.x) max(MESH.x) ...
              min(MESH.y) max(MESH.y) ...
              min(MESH.z) max(MESH.z)]);
        xlabel('x')
        ylabel('y')
        zlabel('z')
        subplot(2,2,3)
        h = slice(MESH.X, MESH.Y, MESH.Z, wf_z, 0, 0, LL.HUB_HT, 'cubic');
        set(h,'FaceColor',          'interp', ...
              'EdgeColor',          'none', ...
              'DiffuseStrength',    0.8)
        colorbar
        hold on
        plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
        axis equal
        axis([min(MESH.x) max(MESH.x) ...
              min(MESH.y) max(MESH.y) ...
              min(MESH.z) max(MESH.z)]);
        xlabel('x')
        ylabel('y')
        zlabel('z')
    subplot(2,2,4)
        h = slice(MESH.X, MESH.Y, MESH.Z, wf_mag, 0, 0, LL.HUB_HT, 'cubic');
        set(h,'FaceColor',          'interp', ...
              'EdgeColor',          'none', ...
              'DiffuseStrength',    0.8)
        colorbar
        hold on
        plot3(xp(1,:), ...
              xp(2,:), ...
              xp(3,:), ...
              '.k');
        axis equal
        axis([min(MESH.x) max(MESH.x) ...
              min(MESH.y) max(MESH.y) ...
              min(MESH.z) max(MESH.z)]);
        xlabel('x')
        ylabel('y')
        zlabel('z')
        drawnow

    %% Vorticity Iso-Surfaces
    figure
    set(gcf,'name', 'Vorticity Isosurfaces','numbertitle','off')
    clf
    isosurface(MESH.X, MESH.Y, MESH.Z, wf_mag, 0.25 * max(max(max(wf_mag))))
    isosurface(MESH.X, MESH.Y, MESH.Z, wf_mag, 0.50 * max(max(max(wf_mag))))
    isosurface(MESH.X, MESH.Y, MESH.Z, wf_mag, 0.75 * max(max(max(wf_mag))))
    hold on
    plot3(xp(1,:), ...
          xp(2,:), ...
          xp(3,:), ...
          '.g');
    alpha(0.25)
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    %% Velocity Field
    figure
    set(gcf,'name', 'Velocity Field','numbertitle','off')
%     subplot(2,2,1)
%     h = slice(MESH.X, MESH.Y, MESH.Z, uf_x, 0, 0, LL.HUB_HT, 'cubic');
%     set(h,'FaceColor',          'interp', ...
%           'EdgeColor',          'none', ...
%           'DiffuseStrength',    0.8)
%     colorbar
%     hold on
%     plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
%     axis equal
%     axis([min(MESH.x) max(MESH.x) ...
%           min(MESH.y) max(MESH.y) ...
%           min(MESH.z) max(MESH.z)]);
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     subplot(2,2,2)
%     h = slice(MESH.X, MESH.Y, MESH.Z, uf_y, 0, 0, LL.HUB_HT, 'cubic');
%     set(h,'FaceColor',          'interp', ...
%           'EdgeColor',          'none', ...
%           'DiffuseStrength',    0.8)
%     colorbar
%     hold on
%     plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
%     axis equal
%     axis([min(MESH.x) max(MESH.x) ...
%           min(MESH.y) max(MESH.y) ...
%           min(MESH.z) max(MESH.z)]);
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     subplot(2,2,3)
%     h = slice(MESH.X, MESH.Y, MESH.Z, uf_z, 0, 0, LL.HUB_HT, 'cubic');
%     set(h,'FaceColor',          'interp', ...
%           'EdgeColor',          'none', ...
%           'DiffuseStrength',    0.8)
%     colorbar
%     hold on
%     plot3(xp(1,:), xp(2,:), xp(3,:), '.k')
%     axis equal
%     axis([min(MESH.x) max(MESH.x) ...
%           min(MESH.y) max(MESH.y) ...
%           min(MESH.z) max(MESH.z)]);
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     subplot(2,2,4)
    h = slice(MESH.X, MESH.Y, MESH.Z, uf_mag, 0, 0, LL.HUB_HT, 'cubic');
    set(h,'FaceColor',          'interp', ...
          'EdgeColor',          'none', ...
          'DiffuseStrength',    0.8)
    colorbar
    hold on
    plot3(xp(1,:), ...
          xp(2,:), ...
          xp(3,:), ...
          '.g');
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')


end

if DEBUG_LVL > 2    
    
    %% Particle Vorticity
    figure
    set(gcf,'name', 'Particle Vorticity','numbertitle','off')
    subplot(2,2,1)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', ap(1,:))
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,2)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', ap(2,:))
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,3)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', ap(3,:))
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    subplot(2,2,4)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', ap_mag)
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')

    % Particle Velocity
    figure
    set(gcf,'name', 'Particle Velocity','numbertitle','off')
    subplot(2,2,1)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(1,:))
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow

    subplot(2,2,2)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(2,:))
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow

    subplot(2,2,3)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(3,:))
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow

    subplot(2,2,4)
    plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up_mag)
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow
end