function [uf, up] = add_freestream(uf, up, ENV)

% velocity field
uf{1} = uf{1} + ENV.velFree(1);
uf{2} = uf{2} + ENV.velFree(2);
uf{3} = uf{3} + ENV.velFree(3);

% particle velocity
up(1,:) = up(1,:) + ENV.velFree(1);
up(2,:) = up(2,:) + ENV.velFree(2);
up(3,:) = up(3,:) + ENV.velFree(3);

end

