% Parallel Arnold tongue map (parfor) for the Winfree-type two-phase model.
% Requires Parallel Computing Toolbox.
%
% x-axis: deltaomega = omega2 - omega1
% y-axis: sigma (0 to 5 by default)
% color : mean angular-velocity ratio (omega2_bar / omega1_bar)
%         computed over final 60 s of each 120 s trial.

clearvars;
close all;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

%% User settings
omega1 = 2.5 * pi;
phi0 = [0, 0.9 * pi];

trial_duration = 120;      % total simulation time [s]
measure_window = 60;       % averaging window from the end [s]
sim_dt = 0.01;              % fixed-step size for deterministic runtime

% Sweep grids
deltaomega_values = linspace(-1, 1, 21);
sigma_values = linspace(0, 5, 21);

% Reconstruction settings (single precompute)
recon_opts = struct();
recon_opts.signal_role = 'a2';
recon_opts.agent_mode = 'phase_id_2';
recon_opts.omega = [omega1, omega1];
recon_opts.sigma = 0;
recon_opts.tspan = [0, 1];
recon_opts.phi0 = phi0;
recon_opts.n_phi_integral = 801;
recon_opts.n_psi_scan = 41;

%% Reconstruct s2 and z_opt from export
recon = simulate_winfree_from_export(recon_opts);
s2_eval = recon.s2_eval;

fprintf('Reconstruction source : %s\n', recon.meta.mat_path);
fprintf('Agent id (s2 fit)     : %d\n', recon.meta.agent_id);
fprintf('psi_plus_opt          : %.12f rad\n', recon.psi_plus_opt);

%% Precompute interpolants for fast evaluation inside ODEs
% Build a coarse grid for s2 and z to avoid expensive basis evaluation in inner loop.
precomp_n = 61; % faster inner loop; increase if you need more accuracy
phi_vals_pre = linspace(0, 2*pi, precomp_n + 1);
phi_vals_pre = phi_vals_pre(1:end-1); % drop duplicate 2*pi endpoint
[Phi1g, Phi2g] = meshgrid(phi_vals_pre, phi_vals_pre);
S2grid = s2_eval(Phi1g, Phi2g);
% gridded interpolant for s2(phi1,phi2)
s2F = griddedInterpolant({phi_vals_pre, phi_vals_pre}, S2grid, 'linear', 'nearest');

% z: sample reconstructed z_opt on same phi grid and build 1D interpolant
z_sample = interp1(recon.z_opt_grid.theta, recon.z_opt_grid.value, phi_vals_pre, 'pchip');
zF = griddedInterpolant(phi_vals_pre', z_sample', 'linear', 'nearest');

fprintf('Using precomputed interpolants: grid %d x %d\n', precomp_n, precomp_n);

%% Arnold map sweep (sequential)
n_sigma = numel(sigma_values);
n_dw = numel(deltaomega_values);
n_cases = n_sigma * n_dw;

progress_every = 5;
fprintf('Starting sequential sweep: %d cases\n', n_cases);
fprintf('Parallel progress: %d/%d (0.00%%), elapsed 0.0s, ETA inf\n', 0, n_cases);
drawnow;

ratio_vec = nan(n_cases, 1);
omega1_bar_vec = nan(n_cases, 1);
omega2_bar_vec = nan(n_cases, 1);

tic;
for i_case = 1:n_cases
    [i_sigma, i_dw] = ind2sub([n_sigma, n_dw], i_case);

    sigma_i = sigma_values(i_sigma);
    dw_i = deltaomega_values(i_dw);
    omega2_i = omega1 + dw_i;

    omega1_bar = NaN;
    omega2_bar = NaN;
    ratio21 = NaN;

    try
        [omega1_bar, omega2_bar, ratio21] = simulate_case_fixed_step(omega1, omega2_i, sigma_i, phi0(:), trial_duration, measure_window, sim_dt, s2F, zF);
    catch ME
        fprintf('Case failed at sigma=%.6g, deltaomega=%.6g: %s\n', sigma_i, dw_i, ME.message);
        drawnow;
    end

    omega1_bar_vec(i_case) = omega1_bar;
    omega2_bar_vec(i_case) = omega2_bar;
    ratio_vec(i_case) = ratio21;

    if mod(i_case, progress_every) == 0 || i_case == 1 || i_case == n_cases
        elapsed_case = toc;
        frac = i_case / n_cases;
        eta = elapsed_case / max(frac, eps) - elapsed_case;
        fprintf('Progress: %d/%d (%.2f%%), elapsed %.1fs, ETA %.1fs\n', i_case, n_cases, 100 * frac, elapsed_case, max(eta, 0));
        drawnow;
    end
end
elapsed = toc;
fprintf('Sequential sweep finished in %.1f s for %d cases.\n', elapsed, n_cases);

ratio_map = reshape(ratio_vec, [n_sigma, n_dw]);
omega1_bar_map = reshape(omega1_bar_vec, [n_sigma, n_dw]);
omega2_bar_map = reshape(omega2_bar_vec, [n_sigma, n_dw]);

%% Plot heatmap
figure('Color', 'w');
imagesc(deltaomega_values, sigma_values, ratio_map);
axis xy;
colormap(turbo);
cb = colorbar;
cb.Label.String = 'mean angular-velocity ratio  \omega_2 / \omega_1  (final 60 s)';
set(gca, 'Layer', 'top');

xlabel('$$\Delta\omega = \omega_2 - \omega_1$$', 'Interpreter', 'latex');
ylabel('$$\sigma$$', 'Interpreter', 'latex');
title('Arnold tongue map (parfor) from reconstructed $$s_2$$ and $$z_{\mathrm{opt}}$$', ...
    'Interpreter', 'latex');

hold on;
lock_tol = 0.01;
lock_mask = isfinite(ratio_map) & abs(ratio_map - 1) <= lock_tol;
h_lock = contour(deltaomega_values, sigma_values, double(lock_mask), [0.5 0.5], 'k-', 'LineWidth', 2.0);
line_handles = findobj(h_lock, 'Type', 'Line');
if numel(line_handles) > 1
    line_lengths = arrayfun(@(h) numel(h.XData), line_handles);
    [~, keep_idx] = max(line_lengths);
    delete(line_handles(setdiff(1:numel(line_handles), keep_idx)));
elseif isempty(line_handles)
    set(h_lock, 'Visible', 'off');
end
if isgraphics(h_lock)
    set(h_lock, 'LineColor', 'k');
end
display_name_set = false;
if isgraphics(h_lock)
    try
        h_lock.DisplayName = '1:1 locking boundary';
        display_name_set = true;
    catch
    end
end
if ~display_name_set
    legend_entry = plot(nan, nan, 'k-', 'LineWidth', 2.0, 'DisplayName', '1:1 locking boundary');
    set(legend_entry, 'Visible', 'off');
end
legend('Location', 'best');

fprintf('Done.\n');

function [omega1_bar, omega2_bar, ratio21] = simulate_case_fixed_step(omega1, omega2, sigma, phi0, trial_duration, window_length, sim_dt, s2F, zF)
    n_steps = max(1, ceil(trial_duration / sim_dt));
    h = trial_duration / n_steps;
    window_start_step = floor((trial_duration - window_length) / h);
    if window_start_step < 0
        window_start_step = 0;
    end

    y = phi0(:);
    y_window_start = y;

    for k = 1:n_steps
        % Forward Euler step: simplest and cheapest time integration.
        f = phase_rhs(y, omega1, omega2, sigma, s2F, zF);
        y = y + h * f;

        if k == window_start_step
            y_window_start = y;
        end
    end

    window_dt = trial_duration - window_start_step * h;
    if window_dt <= 0
        omega1_bar = NaN;
        omega2_bar = NaN;
        ratio21 = NaN;
        return;
    end

    omega1_bar = (y(1) - y_window_start(1)) / window_dt;
    omega2_bar = (y(2) - y_window_start(2)) / window_dt;
    if abs(omega1_bar) < eps
        ratio21 = NaN;
    else
        ratio21 = omega2_bar / omega1_bar;
    end
end

function dy = phase_rhs(y, omega1, omega2, sigma, s2F, zF)
    th1 = mod(y(1), 2*pi);
    th2 = mod(y(2), 2*pi);
    dy = [ ...
        omega1 + sigma * zF(th1) * s2F(th1, th2); ...
        omega2 + sigma * zF(th2) * s2F(th2, th1) ...
        ];
end
