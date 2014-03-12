function output_exampleTurbine( SIM, MESH, PART, ENV, Coord, LL, ...
                                fig1, fig2, ...                                
                                tspan, ...                                
                                xp, ...
                                xp_emit_rotor, ...
                                azim_blades_old, ...
                                azim_blades, ...
                                d_dt_theta_old, ...
                                d_dt_theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here




set(0,'CurrentFigure',fig1);
    subplot(3,1,1)
        hold on
        plot(tspan, [azim_blades_old(1); azim_blades(1)], 'o-r');
        plot(tspan, [azim_blades_old(2); azim_blades(2)], 'o-b');
        plot(tspan, [azim_blades_old(3); azim_blades(3)], 'o-g');
        ylabel('blade azimuth [deg]')
        axis([0 SIM.endtime 0 360]);
        box on
    subplot(3,1,2)
        hold on
        plot(tspan, [    d_dt_theta_old;     d_dt_theta] .* 30/pi, '.-k');
        ylabel('rotor speed [rpm]')
        xlim([0 SIM.endtime]);
        box on
%                 subplot(3,1,3)
%                     hold on
%                     plot(tspan, [    d_dt_theta_old;     d_dt_theta] .* LL.ROTOR_DIA/(2*ENV.velFree), '.-k');
%                     ylabel('tip speed ratio')
%                     xlim([0 SIM.endtime]);
%                     box on
    subplot(3,1,3)
        hold on
        plot(tspan, [PART.nPart_old;        PART.nPart       ] .* LL.ROTOR_DIA/(2*ENV.velFree), 'o-k');
        plot(tspan, [PART.nPart_old_wake;   PART.nPart_wake  ] .* LL.ROTOR_DIA/(2*ENV.velFree), '.-b');
        plot(tspan, [PART.nPart_old_inflow; PART.nPart_inflow] .* LL.ROTOR_DIA/(2*ENV.velFree), '.-r');
        legend('total','wake','inflow');
        ylabel('particles')
        xlabel('time [s]')
        xlim([0 SIM.endtime]);
        box on
drawnow

set(0,'CurrentFigure',fig2);
    clf(fig2);
    hold on
    % plots the global origin (center of rotation) and the
    % energy extraction area (ECA)
    plot3(Coord(1).Origin.Global(1), ...
          Coord(1).Origin.Global(2), ...
          Coord(1).Origin.Global(3), 'ok');                
    % plots the blade lines
    plot3([Coord(1).Origin.Blade(1) Coord(1).Tip.Blade(1)], ...
          [Coord(1).Origin.Blade(2) Coord(1).Tip.Blade(2)], ...
          [Coord(1).Origin.Blade(3) Coord(1).Tip.Blade(3)], 'x-r');
    plot3([Coord(2).Origin.Blade(1) Coord(2).Tip.Blade(1)], ...
          [Coord(2).Origin.Blade(2) Coord(2).Tip.Blade(2)], ...
          [Coord(2).Origin.Blade(3) Coord(2).Tip.Blade(3)], 'x-b');
    plot3([Coord(3).Origin.Blade(1) Coord(3).Tip.Blade(1)], ...
          [Coord(3).Origin.Blade(2) Coord(3).Tip.Blade(2)], ...
          [Coord(3).Origin.Blade(3) Coord(3).Tip.Blade(3)], 'x-g');
    % plots the emitter particles along the blades
    plot3(xp_emit_rotor{1}(:,1), ...
          xp_emit_rotor{1}(:,2), ...
          xp_emit_rotor{1}(:,3), 'or');
    plot3(xp_emit_rotor{2}(:,1), ...
          xp_emit_rotor{2}(:,2), ...
          xp_emit_rotor{2}(:,3), 'ob');
    plot3(xp_emit_rotor{3}(:,1), ...
          xp_emit_rotor{3}(:,2), ...
          xp_emit_rotor{3}(:,3), 'og');
    % plot the particles subject to VTE colored by their strength
    plot3(xp(1,:), ...
          xp(2,:), ...
          xp(3,:), ...
          '.k');
    axis equal
    axis([min(MESH.x) max(MESH.x) ...
          min(MESH.y) max(MESH.y) ...
          min(MESH.z) max(MESH.z)]);
    view(3)
    camproj('perspective')
    box on
    xlabel('streamwise, x [m]')
    ylabel('transverse, y [m]')
    zlabel('depth, z [m]')
drawnow

end % output_exampleTurbine

