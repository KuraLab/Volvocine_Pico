% Arnold tongue map for the Winfree-type two-phase model reconstructed from export.
% x-axis: deltaomega = omega2 - omega1
% y-axis: sigma (0 to 5 by default)
% color : mean angular-velocity ratio (omega2_bar / omega1_bar)
%         computed over final 60 s of each 120 s trial.

clearvars;
close all;

%% User settings
omega1 = 2.5 * pi;
phi0 = [0, 0.9 * pi];

trial_duration = 120;      % total simulation time [s]
measure_window = 60;       % averaging window from the end [s]

% Sweep grids
deltaomega_values = linspace(-1.5 * pi, 1.5 * pi, 101);
sigma_values = linspace(0, 5, 81);

% ODE accuracy (loosen for speed if needed)
ode_rel_tol = 1e-6;
ode_abs_tol = 1e-8;

% Reconstruct settings (run once)
recon_opts = struct();
recon_opts.signal_role = 'a2';
recon_opts.agent_mode = 'phase_id_2';
recon_opts.omega = [omega1, omega1];
recon_opts.sigma = 0;
recon_opts.tspan = [0, 1];
recon_opts.phi0 = phi0;

% These two dominate precompute cost for z_opt extraction.
recon_opts.n_phi_integral = 801;
recon_opts.n_psi_scan = 41;
recon_opts.ode_rel_tol = ode_rel_tol;
recon_opts.ode_abs_tol = ode_abs_tol;

%% Reconstruct s2 and z_opt from export (single precompute)
recon = simulate_winfree_from_export(recon_opts);
s2_eval = recon.s2_eval;
z_eval = recon.z_opt_eval;

fprintf('Reconstruction source : %s\n', recon.meta.mat_path);
fprintf('Agent id (s2 fit)     : %d\n', recon.meta.agent_id);
fprintf('psi_plus_opt          : %.12f rad\n', recon.psi_plus_opt);

%% Arnold map sweep
n_sigma = numel(sigma_values);
n_dw = numel(deltaomega_values);
ratio_map = nan(n_sigma, n_dw);
omega1_bar_map = nan(n_sigma, n_dw);
omega2_bar_map = nan(n_sigma, n_dw);

ode_opts = odeset('RelTol', ode_rel_tol, 'AbsTol', ode_abs_tol);

tic;
for i_sigma = 1:n_sigma
    sigma_i = sigma_values(i_sigma);

    for i_dw = 1:n_dw
        dw_i = deltaomega_values(i_dw);
        omega2_i = omega1 + dw_i;

        odefun = @(t, y) [ ...
            omega1 + sigma_i * z_eval(y(1)) * s2_eval(y(1), y(2)); ...
            omega2_i + sigma_i * z_eval(y(2)) * s2_eval(y(2), y(1)) ...
            ];

        [t, phi] = ode45(odefun, [0, trial_duration], phi0(:), ode_opts);

        [omega1_bar, omega2_bar, ratio21] = compute_final_window_stats(t, phi, measure_window, trial_duration);

        omega1_bar_map(i_sigma, i_dw) = omega1_bar;
        omega2_bar_map(i_sigma, i_dw) = omega2_bar;
        ratio_map(i_sigma, i_dw) = ratio21;
    end

    elapsed = toc;
    frac = i_sigma / n_sigma;
    eta = elapsed / max(frac, eps) - elapsed;
    fprintf('Sweep progress: %3d/%3d (%.1f%%), elapsed %.1fs, ETA %.1fs\n', ...
        i_sigma, n_sigma, 100 * frac, elapsed, max(eta, 0));
end

%% Plot heatmap
figure('Color', 'w');
imagesc(deltaomega_values, sigma_values, ratio_map);
axis xy;
colormap(turbo);
cb = colorbar;
cb.Label.String = 'mean angular-velocity ratio  \omega_2 / \omega_1  (final 60 s)';

xlabel('$$\Delta\omega = \omega_2 - \omega_1$$', 'Interpreter', 'latex');
ylabel('$$\sigma$$', 'Interpreter', 'latex');
title('Arnold tongue map from reconstructed $$s_2$$ and $$z_{\mathrm{opt}}$$', 'Interpreter', 'latex');

% Optional contours help visualizing locking plateaus.
hold on;
contour(deltaomega_values, sigma_values, ratio_map, [1 1], 'k-', 'LineWidth', 1.2, ...
    'DisplayName', '1:1 locking');
legend('Location', 'best');

%% Save outputs
out = struct();
out.meta = struct();
out.meta.generated_at = char(datetime('now'));
out.meta.mat_source = recon.meta.mat_path;
out.meta.agent_id = recon.meta.agent_id;
out.meta.signal_role = recon.meta.signal_role;
out.meta.psi_plus_opt = recon.psi_plus_opt;
out.meta.scale_z_opt = recon.scale_z_opt;
out.meta.omega1 = omega1;
out.meta.phi0 = phi0;
out.meta.trial_duration = trial_duration;
out.meta.measure_window = measure_window;
out.meta.ode_rel_tol = ode_rel_tol;
out.meta.ode_abs_tol = ode_abs_tol;

out.deltaomega_values = deltaomega_values;
out.sigma_values = sigma_values;
out.ratio_map = ratio_map;
out.omega1_bar_map = omega1_bar_map;
out.omega2_bar_map = omega2_bar_map;

assignin('base', 'arnold_tongue_result', out);

fprintf('Done. Result exported to workspace variable: arnold_tongue_result\n');

function [omega1_bar, omega2_bar, ratio21] = compute_final_window_stats(t, phi, window_length, t_end)
start_time = t_end - window_length;
mask = t >= start_time;
idx = find(mask);

if numel(idx) < 2
    % Fallback: use entire trace if solver output is too sparse.
    idx = 1:numel(t);
end

i0 = idx(1);
i1 = idx(end);
dt = t(i1) - t(i0);

if dt <= 0
    omega1_bar = NaN;
    omega2_bar = NaN;
    ratio21 = NaN;
    return;
end

omega1_bar = (phi(i1, 1) - phi(i0, 1)) / dt;
omega2_bar = (phi(i1, 2) - phi(i0, 2)) / dt;
ratio21 = omega2_bar / max(abs(omega1_bar), eps);
end
