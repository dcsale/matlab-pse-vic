%% This function should perform the following tasks:
% 1 - get access to the entire MATLAB workspace during the "main loop", or the output .mat file
% 2 - creates various figures for post-processing or "in-the-loop" sanity check type visualizations.

clear all;
close all;
clc;
load vortex_rings_post_process.mat

% %% inputs - copy and paste
% times1 = [];
% times2 = [];

%% Post Processing


rhs_evals1 = numel(times1);
rhs_evals2 = numel(times2);
rhs_evals3 = numel(times3);
rhs_evals4 = numel(times4);

figure;
set(gcf, 'color', 'white');

hold on; 
plot(times1, 1:rhs_evals1, '-r'); 
plot(times2, 1:rhs_evals2, '-b');
plot(times3, 1:rhs_evals3, '-g');
plot(times4, 1:rhs_evals4, '-c');

xlabel('simulation time, (s)');
ylabel('RHS evaluations');
box on



