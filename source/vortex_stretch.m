function w_du = vortex_stretch(SIM, MESH, wf, uf)
% compute the velocity gradient
[dx_ufx, dy_ufx, dz_ufx] = gradient(uf{1}, MESH.dx(1), MESH.dx(2), MESH.dx(3));
[dx_ufy, dy_ufy, dz_ufy] = gradient(uf{2}, MESH.dx(1), MESH.dx(2), MESH.dx(3));
[dx_ufz, dy_ufz, dz_ufz] = gradient(uf{3}, MESH.dx(1), MESH.dx(2), MESH.dx(3));

% compute the vorticity gradient
% [dx_wfx, dy_wfx, dz_wfx] = gradient(wf{1}, MESH.dx(1), MESH.dx(2), MESH.dx(3));
% [dx_wfy, dy_wfy, dz_wfy] = gradient(wf{2}, MESH.dx(1), MESH.dx(2), MESH.dx(3));
% [dx_wfz, dy_wfz, dz_wfz] = gradient(wf{3}, MESH.dx(1), MESH.dx(2), MESH.dx(3));

% compute the vortex stretching term
w_du_x = wf{1}.*dx_ufx + wf{2}.*dy_ufx + wf{3}.*dz_ufx;
w_du_y = wf{1}.*dx_ufy + wf{2}.*dy_ufy + wf{3}.*dz_ufy;
w_du_z = wf{1}.*dx_ufz + wf{2}.*dy_ufz + wf{3}.*dz_ufz;
w_du   = {w_du_x; w_du_y; w_du_z};

end
