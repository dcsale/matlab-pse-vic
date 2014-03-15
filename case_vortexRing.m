function [ output_args ] = case_vortexRing( input_args )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% =======================================================================%
% update the mesh used for visualization
% ========================================================================%
if MESH.adaptive
    MESH = updateMesh(MESH, RING.Rmajor, xp, PART);
end    

% particle velocity (P2P)
up = vel_P2P(xp, ap, PART, SIM, ENV);

% RHS of the vorticity eqn (P2P)
dap_dt_RHS = get_RHS(PART.nPart, xp, ap, PART.hp, ENV.kin_visc, SIM.runMode_RHS);

% velocity field (P2M)
[uf_x, uf_y, uf_z] = vel_P2M_reg(xp, ap, PART, SIM, MESH, ENV);

% vorticity field (P2M)
[wf_x, wf_y, wf_z] = vort_P2M(MESH, PART.nPart, xp, ap, PART.hp, SIM.runMode_P2M);

% make figures for verification
makeFigures(MESH, RING, xp, wp, up, uf_x, uf_y, uf_z, wf_x, wf_y, wf_z, Inf)

end

