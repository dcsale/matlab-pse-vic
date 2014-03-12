function stop = assignParticleMask(mesh, xp)

stop = false;
if any(xp(1,:) < mesh.x(1)) || any(xp(1,:) > mesh.x(end))
  stop = true;
end
if any(xp(2,:) < mesh.y(1)) || any(xp(2,:) > mesh.y(end))
  stop = true;
end
if any(xp(3,:) < mesh.z(1)) || any(xp(3,:) > mesh.z(end))
  stop = true;
end

end % function