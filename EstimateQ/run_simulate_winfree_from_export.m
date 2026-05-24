% Quick runner for simulate_winfree_from_export
% Run from MATLAB with current folder at EstimateQ (or with EstimateQ on path).

clearvars;
close all;

opts = struct();
opts.signal_role = 'a2';
opts.agent_mode = 'phase_id_2';
opts.omega = [2.5*pi, 2.5*pi];
opts.sigma = 1.00;
opts.tspan = [0, 120];
opts.phi0 = [0, 0.9*pi];

% opts.mat_path can be absolute path or path relative to EstimateQ.
% opts.mat_path = 'Spring2\\omega245\\gamma_exports\\gamma_export_latest.mat';

result = simulate_winfree_from_export(opts);

fprintf('MAT source       : %s\n', result.meta.mat_path);
fprintf('Agent id (s2 fit): %d\n', result.meta.agent_id);
fprintf('psi_plus_opt     : %.12f rad\n', result.psi_plus_opt);
fprintf('scale_z_opt      : %.12f\n', result.scale_z_opt);

phi_wrapped = mod(result.phi, 2 * pi);
phase_diff = mod(phi_wrapped(:, 1) - phi_wrapped(:, 2) + pi, 2 * pi) - pi;

figure('Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
plot(result.t, phi_wrapped(:, 1), 'LineWidth', 1.4, 'DisplayName', '$$\phi_1\ \mathrm{mod}\ 2\pi$$');
hold on;
plot(result.t, phi_wrapped(:, 2), 'LineWidth', 1.4, 'DisplayName', '$$\phi_2\ \mathrm{mod}\ 2\pi$$');
grid on;
box on;
ylabel('phase [rad]');
legend('Location', 'best', 'Interpreter', 'latex');
title('Winfree model simulation from exported s_2 and z_{opt}');

nexttile;
plot(result.t, phase_diff, 'LineWidth', 1.6);
grid on;
box on;
ylabel('$$\phi_1-\phi_2$$', 'Interpreter', 'latex');

nexttile;
plot(result.z_opt_grid.theta, result.z_opt_grid.value, 'LineWidth', 1.6);
grid on;
box on;
xlabel('$$\theta$$', 'Interpreter', 'latex');
ylabel('$$z_{\mathrm{opt}}(\theta)$$', 'Interpreter', 'latex');
xticks([0, pi/2, pi, 3*pi/2, 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});

% Optional: quick snapshot of reconstructed s2(phi1, phi2)
n_grid = 121;
phi_vals = linspace(0, 2 * pi, n_grid);
[Phi1, Phi2] = meshgrid(phi_vals, phi_vals);
S2 = result.s2_eval(Phi1, Phi2);

figure('Color', 'w');
imagesc(phi_vals, phi_vals, S2);
axis xy;
axis equal;
colormap(turbo);
colorbar;
xlabel('$$\phi_1$$', 'Interpreter', 'latex');
ylabel('$$\phi_2$$', 'Interpreter', 'latex');
title('Reconstructed $$s_2(\phi_1,\phi_2)$$', 'Interpreter', 'latex');
xticks([0, pi/2, pi, 3*pi/2, 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
yticks([0, pi/2, pi, 3*pi/2, 2*pi]);
yticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
