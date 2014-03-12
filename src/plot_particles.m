function plot_particles(xp, wp, MESH, tag)

wp_mag = sqrt( wp(1,:).^2 + wp(2,:).^2 + wp(3,:).^2 );

%% Particle quantity
figure
set(gcf,'name', tag, 'numbertitle', 'off') 

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
plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', wp_mag)
axis equal
axis([MESH.xmin(1) MESH.xmax(1) ...
      MESH.xmin(2) MESH.xmax(2) ...
      MESH.xmin(3) MESH.xmax(3)]);
xlabel('x')
ylabel('y')
zlabel('z')

% Particle Velocity
% figure
% set(gcf,'name', 'Particle Velocity','numbertitle','off')
% subplot(2,2,1)
% plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(1,:))
% axis equal
% axis([min(MESH.x) max(MESH.x) ...
%       min(MESH.y) max(MESH.y) ...
%       min(MESH.z) max(MESH.z)]);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% drawnow
% 
% subplot(2,2,2)
% plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(2,:))
% axis equal
% axis([min(MESH.x) max(MESH.x) ...
%       min(MESH.y) max(MESH.y) ...
%       min(MESH.z) max(MESH.z)]);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% drawnow
% 
% subplot(2,2,3)
% plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up(3,:))
% axis equal
% axis([min(MESH.x) max(MESH.x) ...
%       min(MESH.y) max(MESH.y) ...
%       min(MESH.z) max(MESH.z)]);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% drawnow
% 
% subplot(2,2,4)
% plot3k({xp(1,:) xp(2,:) xp(3,:)}, 'ColorData', up_mag)
% axis equal
% axis([min(MESH.x) max(MESH.x) ...
%       min(MESH.y) max(MESH.y) ...
%       min(MESH.z) max(MESH.z)]);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% drawnow

end % plot_particles

