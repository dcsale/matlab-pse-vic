%% =======================================================================%
% benchmarking Poisson solver
% ========================================================================%
clc
NX = [64 128 256 288 320]; % 384 is too big and uses 64GB + swap

for n = 1:numel(NX)
   MESH.NX = [NX(n), NX(n), NX(n)];

   MESH = init_mesh(SIM, MESH);                   % init mesh
   wf   = init_field(SIM, MESH);                  % init vorticity field
   uf   = PoissonSolve3D(wf, MESH, SIM, true);          % solve Poisson eqn for velocity

   clear MESH wf uf

end