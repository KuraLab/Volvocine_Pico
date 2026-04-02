% Minimal script: reconstruct and plot fitted Fourier surface
% from gamma_export_latest.mat.

clearvars;
close all;

% --- User settings ---
mat_path = '';          % '' -> auto-pick newest gamma_export_latest.mat under EstimateQ
signal_role = 'a2';     % 'a2' or 'derived'
n_grid = 121;
n_phi_line = 401;       % sampling count for 1D W(phi) plot
n_phi_integral = 2001;  % sampling count for numerical integral over phi
n_psi_scan = 101;       % sampling count for psi_+ sweep in R(psi_+)
psi_scan_min = -pi;
psi_scan_max = pi;
export_fourier_params = true;  % true -> print/save fitted Fourier parameters
fourier_param_output_dir = ''; % '' -> use the folder that contains mat_path
processing_modes = struct('s1', false, 's2', true);  % enable/disable per-s processing

if ~processing_modes.s1 && ~processing_modes.s2
    error('At least one of processing_modes.s1 or processing_modes.s2 must be true.');
end

if isempty(mat_path)
    files = dir(fullfile('EstimateQ', 'Spring*', '*', 'gamma_exports', 'gamma_export_latest.mat'));
    assert(~isempty(files), 'No gamma_export_latest.mat was found under EstimateQ.');
    [~, newest_idx] = max([files.datenum]);
    mat_path = fullfile(files(newest_idx).folder, files(newest_idx).name);
end

S = load(mat_path, 'gamma_export');
assert(isfield(S, 'gamma_export'), 'The MAT file does not contain gamma_export.');
G = S.gamma_export;

assert(isfield(G, 'phase_agent_ids') && numel(G.phase_agent_ids) >= 2, ...
    'gamma_export.phase_agent_ids does not contain two phase-agent IDs.');
phase_ids = G.phase_agent_ids(:).';
agent_id_s1 = phase_ids(1);
agent_id_s2 = phase_ids(2);

fit_data_s1 = get_fit_data_for_agent(G, agent_id_s1, signal_role);
fit_data_s2 = get_fit_data_for_agent(G, agent_id_s2, signal_role);

coeff_s1 = [];
basis_types_s1 = {};
m_order_s1 = [];
n_order_s1 = [];
z_mean_s1 = NaN;
if processing_modes.s1
    [coeff_s1, basis_types_s1, m_order_s1, n_order_s1, z_mean_s1] = extract_basis_data(fit_data_s1);
end

coeff_s2 = [];
basis_types_s2 = {};
m_order_s2 = [];
n_order_s2 = [];
z_mean_s2 = NaN;
if processing_modes.s2
    [coeff_s2, basis_types_s2, m_order_s2, n_order_s2, z_mean_s2] = extract_basis_data(fit_data_s2);
end

phi_values = linspace(0, 2*pi, n_grid);
[Phi1, Phi2] = meshgrid(phi_values, phi_values);

if processing_modes.s1
    s1_values = evaluate_exported_s(Phi1, Phi2, coeff_s1, basis_types_s1, m_order_s1, n_order_s1) + z_mean_s1;

    figure('Color', 'w');
    surf(Phi1, Phi2, s1_values, 'EdgeColor', 'none');
    view(40, 30);
    grid on;
    box on;
    axis tight;
    colormap(turbo);
    colorbar;
    xlabel('\phi_1');
    ylabel('\phi_2');
    zlabel('s(\phi_1,\phi_2)');
    title(sprintf('s_1(\\phi_1,\\phi_2) | agent %d | %s', agent_id_s1, mat_path), 'Interpreter', 'none');
    tuneFigure;
end

if processing_modes.s2
    s2_values = evaluate_exported_s(Phi1, Phi2, coeff_s2, basis_types_s2, m_order_s2, n_order_s2) + z_mean_s2;

    figure('Color', 'w');
    surf(Phi1, Phi2, s2_values, 'EdgeColor', 'none');
    view(40, 30);
    grid on;
    box on;
    axis tight;
    colormap(turbo);
    colorbar;
    xlabel('\phi_1');
    ylabel('\phi_2');
    zlabel('s(\phi_1,\phi_2)');
    title(sprintf('s_2(\\phi_1,\\phi_2) | agent %d | %s', agent_id_s2, mat_path), 'Interpreter', 'none');
    tuneFigure;
end

% R(psi_+) = integral_0^{2*pi} W(phi)^2 dphi
phi_integral = linspace(0, 2*pi, n_phi_integral);
psi_scan_values = linspace(psi_scan_min, psi_scan_max, n_psi_scan);
R1_values = nan(size(psi_scan_values));
R2_values = nan(size(psi_scan_values));

for i = 1:numel(psi_scan_values)
    psi_i = psi_scan_values(i);

    if processing_modes.s1
        W1_i = compute_w_profile(phi_integral, psi_i, coeff_s1, basis_types_s1, m_order_s1, n_order_s1, 's1');
        R1_values(i) = trapz(phi_integral, W1_i .^ 2);
    end

    if processing_modes.s2
        W2_i = compute_w_profile(phi_integral, psi_i, coeff_s2, basis_types_s2, m_order_s2, n_order_s2, 's2');
        R2_values(i) = trapz(phi_integral, W2_i .^ 2);
    end
end

R1_max = NaN;
R2_max = NaN;
idx_max_s1 = NaN;
idx_max_s2 = NaN;
idx_max_s1_alt = NaN;
idx_max_s2_alt = NaN;
psi_plus_s1 = NaN;
psi_plus_s2 = NaN;
psi_plus_s1_alt = NaN;
psi_plus_s2_alt = NaN;
has_symmetric_peak_s1 = false;
has_symmetric_peak_s2 = false;
R1_alt = NaN;
R2_alt = NaN;
if processing_modes.s1
    [R1_max, idx_max_s1, idx_max_s1_alt, has_symmetric_peak_s1, R1_alt] = find_symmetric_max_candidate(R1_values, psi_scan_values);
    psi_plus_s1 = psi_scan_values(idx_max_s1);
    if has_symmetric_peak_s1
        psi_plus_s1_alt = psi_scan_values(idx_max_s1_alt);
    end
end
if processing_modes.s2
    [R2_max, idx_max_s2, idx_max_s2_alt, has_symmetric_peak_s2, R2_alt] = find_symmetric_max_candidate(R2_values, psi_scan_values);
    psi_plus_s2 = psi_scan_values(idx_max_s2);
    if has_symmetric_peak_s2
        psi_plus_s2_alt = psi_scan_values(idx_max_s2_alt);
    end
end

% s1: W(phi) = s(phi, phi-psi_+) - s(phi, phi+psi_+)
% s2: W(phi) = s(phi-psi_+, phi) - s(phi+psi_+, phi)
phi_line = linspace(-pi, pi, n_phi_line);

w1_values = [];
w1_values_alt = [];
if processing_modes.s1
    w1_values = compute_w_profile(phi_line, psi_plus_s1, coeff_s1, basis_types_s1, m_order_s1, n_order_s1, 's1');
    if has_symmetric_peak_s1
        w1_values_alt = compute_w_profile(phi_line, psi_plus_s1_alt, coeff_s1, basis_types_s1, m_order_s1, n_order_s1, 's1');
    end
end

w2_values = [];
w2_values_alt = [];
if processing_modes.s2
    w2_values = compute_w_profile(phi_line, psi_plus_s2, coeff_s2, basis_types_s2, m_order_s2, n_order_s2, 's2');
    if has_symmetric_peak_s2
        w2_values_alt = compute_w_profile(phi_line, psi_plus_s2_alt, coeff_s2, basis_types_s2, m_order_s2, n_order_s2, 's2');
    end
end

% Normalize W power to match unit sine-wave power on the same phi grid.
target_power = trapz(phi_line, sin(phi_line) .^ 2) / (2 * pi);
w1_power_raw = NaN;
w2_power_raw = NaN;
w1_power_raw_alt = NaN;
w2_power_raw_alt = NaN;
if processing_modes.s1
    w1_power_raw = trapz(phi_line, w1_values .^ 2) / (2 * pi);
    if has_symmetric_peak_s1
        w1_power_raw_alt = trapz(phi_line, w1_values_alt .^ 2) / (2 * pi);
    end
end
if processing_modes.s2
    w2_power_raw = trapz(phi_line, w2_values .^ 2) / (2 * pi);
    if has_symmetric_peak_s2
        w2_power_raw_alt = trapz(phi_line, w2_values_alt .^ 2) / (2 * pi);
    end
end

scale_w1 = NaN;
scale_w2 = NaN;
scale_w1_alt = NaN;
scale_w2_alt = NaN;
if processing_modes.s1
    scale_w1 = sqrt(target_power / max(w1_power_raw, eps));
    w1_values = scale_w1 * w1_values;
    if has_symmetric_peak_s1
        scale_w1_alt = sqrt(target_power / max(w1_power_raw_alt, eps));
        w1_values_alt = scale_w1_alt * w1_values_alt;
    end
end
if processing_modes.s2
    scale_w2 = sqrt(target_power / max(w2_power_raw, eps));
    w2_values = scale_w2 * w2_values;
    if has_symmetric_peak_s2
        scale_w2_alt = sqrt(target_power / max(w2_power_raw_alt, eps));
        w2_values_alt = scale_w2_alt * w2_values_alt;
    end
end

w1_power_norm = NaN;
w2_power_norm = NaN;
w1_power_norm_alt = NaN;
w2_power_norm_alt = NaN;
if processing_modes.s1
    w1_power_norm = trapz(phi_line, w1_values .^ 2) / (2 * pi);
    if has_symmetric_peak_s1
        w1_power_norm_alt = trapz(phi_line, w1_values_alt .^ 2) / (2 * pi);
    end
end
if processing_modes.s2
    w2_power_norm = trapz(phi_line, w2_values .^ 2) / (2 * pi);
    if has_symmetric_peak_s2
        w2_power_norm_alt = trapz(phi_line, w2_values_alt .^ 2) / (2 * pi);
    end
end

% Fit W(phi) in Figure 3 by Fourier series up to original s-function order.
fit_order_s1 = NaN;
fit_order_s2 = NaN;
fourier_fit_w1 = struct();
fourier_fit_w2 = struct();
fourier_fit_w1_alt = struct();
fourier_fit_w2_alt = struct();
fourier_param_table_w1 = table();
fourier_param_table_w2 = table();
fourier_param_table_w1_alt = table();
fourier_param_table_w2_alt = table();
if processing_modes.s1
    fit_order_s1 = max(0, round(max([m_order_s1(:); n_order_s1(:)])));
    fourier_fit_w1 = fit_fourier_series_periodic(phi_line, w1_values, fit_order_s1);
    fourier_param_table_w1 = build_fourier_param_table(fourier_fit_w1);
    if has_symmetric_peak_s1
        fourier_fit_w1_alt = fit_fourier_series_periodic(phi_line, w1_values_alt, fit_order_s1);
        fourier_param_table_w1_alt = build_fourier_param_table(fourier_fit_w1_alt);
    end
end
if processing_modes.s2
    fit_order_s2 = max(0, round(max([m_order_s2(:); n_order_s2(:)])));
    fourier_fit_w2 = fit_fourier_series_periodic(phi_line, w2_values, fit_order_s2);
    fourier_param_table_w2 = build_fourier_param_table(fourier_fit_w2);
    if has_symmetric_peak_s2
        fourier_fit_w2_alt = fit_fourier_series_periodic(phi_line, w2_values_alt, fit_order_s2);
        fourier_param_table_w2_alt = build_fourier_param_table(fourier_fit_w2_alt);
    end
end

if export_fourier_params
    if isempty(fourier_param_output_dir)
        output_dir = fileparts(mat_path);
    else
        output_dir = fourier_param_output_dir;
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    mat_path_fit = fullfile(output_dir, sprintf('W_fourier_fit_params_%s.mat', signal_role));

    if processing_modes.s1
        csv_path_w1 = fullfile(output_dir, sprintf('W1_fourier_fit_params_agent%d_%s.csv', agent_id_s1, signal_role));
        writetable(fourier_param_table_w1, csv_path_w1);
    end
    if processing_modes.s2
        csv_path_w2 = fullfile(output_dir, sprintf('W2_fourier_fit_params_agent%d_%s.csv', agent_id_s2, signal_role));
        writetable(fourier_param_table_w2, csv_path_w2);
    end

    fourier_fit_export = struct();
    fourier_fit_export.mat_source = mat_path;
    fourier_fit_export.signal_role = signal_role;
    fourier_fit_export.agent_id_s1 = agent_id_s1;
    fourier_fit_export.agent_id_s2 = agent_id_s2;
    if processing_modes.s1
        fourier_fit_export.psi_plus_s1 = psi_plus_s1;
        fourier_fit_export.fit_order_s1 = fit_order_s1;
        fourier_fit_export.fit_w1 = fourier_fit_w1;
        fourier_fit_export.table_w1 = fourier_param_table_w1;
    end
    if processing_modes.s2
        fourier_fit_export.psi_plus_s2 = psi_plus_s2;
        fourier_fit_export.fit_order_s2 = fit_order_s2;
        fourier_fit_export.fit_w2 = fourier_fit_w2;
        fourier_fit_export.table_w2 = fourier_param_table_w2;
    end
    save(mat_path_fit, 'fourier_fit_export');
end

reference_cos = cos(phi_line + pi - 0.6 * pi);

% z(phi) candidates for Gamma(psi):
% 1) power-normalized W at adopted psi_+^*
% 2) reference cosine used in comparison plot
z1_candidate_from_w = [];
z2_candidate_from_w = [];
z1_candidate_from_w_alt = [];
z2_candidate_from_w_alt = [];
if processing_modes.s1
    z1_candidate_from_w = scale_w1 * compute_w_profile(phi_integral, psi_plus_s1, coeff_s1, basis_types_s1, m_order_s1, n_order_s1, 's1');
    if has_symmetric_peak_s1
        z1_candidate_from_w_alt = scale_w1_alt * compute_w_profile(phi_integral, psi_plus_s1_alt, coeff_s1, basis_types_s1, m_order_s1, n_order_s1, 's1');
    end
end
if processing_modes.s2
    z2_candidate_from_w = scale_w2 * compute_w_profile(phi_integral, psi_plus_s2, coeff_s2, basis_types_s2, m_order_s2, n_order_s2, 's2');
    if has_symmetric_peak_s2
        z2_candidate_from_w_alt = scale_w2_alt * compute_w_profile(phi_integral, psi_plus_s2_alt, coeff_s2, basis_types_s2, m_order_s2, n_order_s2, 's2');
    end
end
z_candidate_sine = cos(phi_integral + pi - 0.6 * pi);

Gamma1_from_w = nan(size(psi_scan_values));
Gamma1_from_sine = nan(size(psi_scan_values));
Gamma2_from_w = nan(size(psi_scan_values));
Gamma2_from_sine = nan(size(psi_scan_values));
Gamma1_from_w_alt = nan(size(psi_scan_values));
Gamma2_from_w_alt = nan(size(psi_scan_values));
if processing_modes.s1
    Gamma1_from_w = compute_gamma_curve_from_z(psi_scan_values, phi_integral, z1_candidate_from_w, ...
        coeff_s1, basis_types_s1, m_order_s1, n_order_s1, z_mean_s1, 's1');
    Gamma1_from_sine = compute_gamma_curve_from_z(psi_scan_values, phi_integral, z_candidate_sine, ...
        coeff_s1, basis_types_s1, m_order_s1, n_order_s1, z_mean_s1, 's1');
    if has_symmetric_peak_s1
        Gamma1_from_w_alt = compute_gamma_curve_from_z(psi_scan_values, phi_integral, z1_candidate_from_w_alt, ...
            coeff_s1, basis_types_s1, m_order_s1, n_order_s1, z_mean_s1, 's1');
    end
end
if processing_modes.s2
    Gamma2_from_w = compute_gamma_curve_from_z(psi_scan_values, phi_integral, z2_candidate_from_w, ...
        coeff_s2, basis_types_s2, m_order_s2, n_order_s2, z_mean_s2, 's2');
    Gamma2_from_sine = compute_gamma_curve_from_z(psi_scan_values, phi_integral, z_candidate_sine, ...
        coeff_s2, basis_types_s2, m_order_s2, n_order_s2, z_mean_s2, 's2');
    if has_symmetric_peak_s2
        Gamma2_from_w_alt = compute_gamma_curve_from_z(psi_scan_values, phi_integral, z2_candidate_from_w_alt, ...
            coeff_s2, basis_types_s2, m_order_s2, n_order_s2, z_mean_s2, 's2');
    end
end

% Odd component: Gamma(psi) - Gamma(-psi)
if psi_scan_min >= 0 || psi_scan_max <= 0
    warning(['psi scan range does not include both signs. ', ...
        'Gamma(-psi) is evaluated by extrapolation and may be inaccurate.']);
end

Gamma1_from_w_neg = nan(size(psi_scan_values));
Gamma1_from_sine_neg = nan(size(psi_scan_values));
Gamma2_from_w_neg = nan(size(psi_scan_values));
Gamma2_from_sine_neg = nan(size(psi_scan_values));
Gamma1_from_w_alt_neg = nan(size(psi_scan_values));
Gamma2_from_w_alt_neg = nan(size(psi_scan_values));

if processing_modes.s1
    Gamma1_from_w_neg = interp1(psi_scan_values, Gamma1_from_w, -psi_scan_values, 'pchip', 'extrap');
    Gamma1_from_sine_neg = interp1(psi_scan_values, Gamma1_from_sine, -psi_scan_values, 'pchip', 'extrap');
    if has_symmetric_peak_s1
        Gamma1_from_w_alt_neg = interp1(psi_scan_values, Gamma1_from_w_alt, -psi_scan_values, 'pchip', 'extrap');
    end
end
if processing_modes.s2
    Gamma2_from_w_neg = interp1(psi_scan_values, Gamma2_from_w, -psi_scan_values, 'pchip', 'extrap');
    Gamma2_from_sine_neg = interp1(psi_scan_values, Gamma2_from_sine, -psi_scan_values, 'pchip', 'extrap');
    if has_symmetric_peak_s2
        Gamma2_from_w_alt_neg = interp1(psi_scan_values, Gamma2_from_w_alt, -psi_scan_values, 'pchip', 'extrap');
    end
end

Gamma1_odd_from_w = Gamma1_from_w - Gamma1_from_w_neg;
Gamma1_odd_from_sine = Gamma1_from_sine - Gamma1_from_sine_neg;
Gamma2_odd_from_w = Gamma2_from_w - Gamma2_from_w_neg;
Gamma2_odd_from_sine = Gamma2_from_sine - Gamma2_from_sine_neg;
Gamma1_odd_from_w_alt = Gamma1_from_w_alt - Gamma1_from_w_alt_neg;
Gamma2_odd_from_w_alt = Gamma2_from_w_alt - Gamma2_from_w_alt_neg;

n_w_tiles = double(processing_modes.s1) + double(processing_modes.s2);
if n_w_tiles > 0
    figure('Color', 'w');
    tiledlayout(n_w_tiles, 1, 'TileSpacing', 'compact');

    if processing_modes.s1
        nexttile;
        plot(phi_line, w1_values, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('W_1, \\psi_+=%.6g', psi_plus_s1));
        hold on;
        plot(phi_line, fourier_fit_w1.y_fit, '-.', 'LineWidth', 1.6, ...
            'DisplayName', sprintf('Fourier fit W_1 (N=%d, \\psi_+=%.6g)', fit_order_s1, psi_plus_s1));
        if has_symmetric_peak_s1
            plot(phi_line, w1_values_alt, 'LineWidth', 1.4, ...
                'DisplayName', sprintf('W_1, \\psi_+=%.6g', psi_plus_s1_alt));
            plot(phi_line, fourier_fit_w1_alt.y_fit, '-.', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Fourier fit W_1 (N=%d, \\psi_+=%.6g)', fit_order_s1, psi_plus_s1_alt));
        end
        plot(phi_line, reference_cos, '--', 'LineWidth', 1.6, 'DisplayName', 'cos(\phi+\pi-0.6\pi)');
        grid on;
        box on;
        xlim([-pi, pi]);
        xticks([-pi, -pi/2, 0, pi/2, pi]);
        xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        xlabel('\phi');
        ylabel('value');
        legend('Location', 'best');
        title(sprintf('s_1: W_1(\\phi) vs cos(\\phi+\\pi-0.6\\pi), \\psi_+^* = %.6g rad', psi_plus_s1));
    end

    if processing_modes.s2
        nexttile;
        plot(phi_line, w2_values, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('W_2, \\psi_+=%.6g', psi_plus_s2));
        hold on;
        plot(phi_line, fourier_fit_w2.y_fit, '-.', 'LineWidth', 1.6, ...
            'DisplayName', sprintf('Fourier fit W_2 (N=%d, \\psi_+=%.6g)', fit_order_s2, psi_plus_s2));
        if has_symmetric_peak_s2
            plot(phi_line, w2_values_alt, 'LineWidth', 1.4, ...
                'DisplayName', sprintf('W_2, \\psi_+=%.6g', psi_plus_s2_alt));
            plot(phi_line, fourier_fit_w2_alt.y_fit, '-.', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Fourier fit W_2 (N=%d, \\psi_+=%.6g)', fit_order_s2, psi_plus_s2_alt));
        end
        plot(phi_line, reference_cos, '--', 'LineWidth', 1.6, 'DisplayName', 'cos(\phi+\pi-0.6\pi)');
        grid on;
        box on;
        xlim([-pi, pi]);
        xticks([-pi, -pi/2, 0, pi/2, pi]);
        xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        xlabel('\phi');
        ylabel('value');
        legend('Location', 'best');
        title(sprintf('s_2: W_2(\\phi) vs cos(\\phi+\\pi-0.6\\pi), \\psi_+^* = %.6g rad', psi_plus_s2));
    end
    tuneFigure;
end

figure('Color', 'w');
hold on;
if processing_modes.s1
    plot(psi_scan_values, R1_values, 'LineWidth', 1.8, 'DisplayName', 'R_1(\psi_+) from W_1');
    plot(psi_plus_s1, R1_max, 'o', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'R_1 max');
    if has_symmetric_peak_s1
        plot(psi_plus_s1_alt, R1_alt, 's', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'R_1 symmetric max');
    end
end
if processing_modes.s2
    plot(psi_scan_values, R2_values, 'LineWidth', 1.8, 'DisplayName', 'R_2(\psi_+) from W_2');
    plot(psi_plus_s2, R2_max, 'o', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'R_2 max');
    if has_symmetric_peak_s2
        plot(psi_plus_s2_alt, R2_alt, 's', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'R_2 symmetric max');
    end
end
grid on;
box on;
xlim([psi_scan_min, psi_scan_max]);
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xlabel('\psi_+');
ylabel('R(\psi_+)');
legend('Location', 'best');
title('R(\psi_+) = \int_0^{2\pi} W(\phi)^2 d\phi');
tuneFigure;

if processing_modes.s1
    figure('Color', 'w');
    plot(psi_scan_values, Gamma1_from_w, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('z_1(\\phi)=W_1(\\phi;\\psi_+=%.6g)', psi_plus_s1));
    hold on;
    if has_symmetric_peak_s1
        plot(psi_scan_values, Gamma1_from_w_alt, 'LineWidth', 1.4, ...
            'DisplayName', sprintf('z_1(\\phi)=W_1(\\phi;\\psi_+=%.6g)', psi_plus_s1_alt));
    end
    plot(psi_scan_values, Gamma1_from_sine, '--', 'LineWidth', 1.8, ...
        'DisplayName', 'z(\phi)=cos(\phi+\pi-0.6\pi)');
    grid on;
    box on;
    xlim([psi_scan_min, psi_scan_max]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    xlabel('\psi');
    ylabel('\Gamma_1(\psi)');
    legend('Location', 'best');
    title('\Gamma_1(\psi) = (1/(2\pi))\int_0^{2\pi} z(\phi)s_1(\phi,\phi-\psi)d\phi');
    tuneFigure;
end

if processing_modes.s2
    figure('Color', 'w');
    plot(psi_scan_values, Gamma2_from_w, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('z_2(\\phi)=W_2(\\phi;\\psi_+=%.6g)', psi_plus_s2));
    hold on;
    if has_symmetric_peak_s2
        plot(psi_scan_values, Gamma2_from_w_alt, 'LineWidth', 1.4, ...
            'DisplayName', sprintf('z_2(\\phi)=W_2(\\phi;\\psi_+=%.6g)', psi_plus_s2_alt));
    end
    plot(psi_scan_values, Gamma2_from_sine, '--', 'LineWidth', 1.8, ...
        'DisplayName', 'z(\phi)=cos(\phi+\pi-0.6\pi)');
    grid on;
    box on;
    xlim([psi_scan_min, psi_scan_max]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    xlabel('\psi');
    ylabel('\Gamma_2(\psi)');
    legend('Location', 'best');
    title('\Gamma_2(\psi) = (1/(2\pi))\int_0^{2\pi} z(\phi)s_2(\phi-\psi,\phi)d\phi');
    tuneFigure;
end

if processing_modes.s1
    figure('Color', 'w');
    plot(psi_scan_values, Gamma1_odd_from_w, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('\\Gamma_1(\\psi)-\\Gamma_1(-\\psi), z_1=W_1(\\psi_+=%.6g)', psi_plus_s1));
    hold on;
    if has_symmetric_peak_s1
        plot(psi_scan_values, Gamma1_odd_from_w_alt, 'LineWidth', 1.4, ...
            'DisplayName', sprintf('\\Gamma_1(\\psi)-\\Gamma_1(-\\psi), z_1=W_1(\\psi_+=%.6g)', psi_plus_s1_alt));
    end
    plot(psi_scan_values, Gamma1_odd_from_sine, '--', 'LineWidth', 1.8, ...
        'DisplayName', '\Gamma_1(\psi)-\Gamma_1(-\psi), z=cos');
    grid on;
    box on;
    xlim([psi_scan_min, psi_scan_max]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    xlabel('\psi');
    ylabel('\Delta\Gamma_1(\psi)');
    legend('Location', 'best');
    title('\Delta\Gamma_1(\psi)=\Gamma_1(\psi)-\Gamma_1(-\psi)');
    tuneFigure;
end

if processing_modes.s2
    figure('Color', 'w');
    plot(psi_scan_values, Gamma2_odd_from_w, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('\\Gamma_2(\\psi)-\\Gamma_2(-\\psi), z_2=W_2(\\psi_+=%.6g)', psi_plus_s2));
    hold on;
    if has_symmetric_peak_s2
        plot(psi_scan_values, Gamma2_odd_from_w_alt, 'LineWidth', 1.4, ...
            'DisplayName', sprintf('\\Gamma_2(\\psi)-\\Gamma_2(-\\psi), z_2=W_2(\\psi_+=%.6g)', psi_plus_s2_alt));
    end
    plot(psi_scan_values, Gamma2_odd_from_sine, '--', 'LineWidth', 1.8, ...
        'DisplayName', '\Gamma_2(\psi)-\Gamma_2(-\psi), z=cos');
    grid on;
    box on;
    xlim([psi_scan_min, psi_scan_max]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    xlabel('\psi');
    ylabel('\Delta\Gamma_2(\psi)');
    legend('Location', 'best');
    title('\Delta\Gamma_2(\psi)=\Gamma_2(\psi)-\Gamma_2(-\psi)');
    tuneFigure;
end

if processing_modes.s1
    fprintf('R_1 max = %.12f at psi_plus^* = %.12f rad\n', R1_max, psi_plus_s1);
    if has_symmetric_peak_s1
        fprintf('R_1 symmetric max = %.12f at psi_plus^* = %.12f rad (sum with primary: %.12e)\n', ...
            R1_alt, psi_plus_s1_alt, psi_plus_s1 + psi_plus_s1_alt);
    else
        fprintf('R_1 symmetric max near -psi_plus^* was not found on the current scan grid.\n');
    end
end
if processing_modes.s2
    fprintf('R_2 max = %.12f at psi_plus^* = %.12f rad\n', R2_max, psi_plus_s2);
    if has_symmetric_peak_s2
        fprintf('R_2 symmetric max = %.12f at psi_plus^* = %.12f rad (sum with primary: %.12e)\n', ...
            R2_alt, psi_plus_s2_alt, psi_plus_s2 + psi_plus_s2_alt);
    else
        fprintf('R_2 symmetric max near -psi_plus^* was not found on the current scan grid.\n');
    end
end
if processing_modes.s1
    fprintf('s_1 baseline z_mean (added back) = %.12f\n', z_mean_s1);
end
if processing_modes.s2
    fprintf('s_2 baseline z_mean (added back) = %.12f\n', z_mean_s2);
end
if processing_modes.s1
    fprintf('W_1 raw power = %.12f, normalized power = %.12f (target sine power = %.12f)\n', ...
        w1_power_raw, w1_power_norm, target_power);
    if has_symmetric_peak_s1
        fprintf('W_1(raw/normalized) for symmetric psi = %.12f / %.12f\n', ...
            w1_power_raw_alt, w1_power_norm_alt);
    end
end
if processing_modes.s2
    fprintf('W_2 raw power = %.12f, normalized power = %.12f (target sine power = %.12f)\n', ...
        w2_power_raw, w2_power_norm, target_power);
    if has_symmetric_peak_s2
        fprintf('W_2(raw/normalized) for symmetric psi = %.12f / %.12f\n', ...
            w2_power_raw_alt, w2_power_norm_alt);
    end
end
if processing_modes.s1
    fprintf('Fourier fit W_1 order = %d, RMSE = %.12f, R^2 = %.12f\n', ...
        fourier_fit_w1.order, fourier_fit_w1.rmse, fourier_fit_w1.r2);
    if has_symmetric_peak_s1
        fprintf('Fourier fit W_1 (symmetric psi) order = %d, RMSE = %.12f, R^2 = %.12f\n', ...
            fourier_fit_w1_alt.order, fourier_fit_w1_alt.rmse, fourier_fit_w1_alt.r2);
    end
end
if processing_modes.s2
    fprintf('Fourier fit W_2 order = %d, RMSE = %.12f, R^2 = %.12f\n', ...
        fourier_fit_w2.order, fourier_fit_w2.rmse, fourier_fit_w2.r2);
    if has_symmetric_peak_s2
        fprintf('Fourier fit W_2 (symmetric psi) order = %d, RMSE = %.12f, R^2 = %.12f\n', ...
            fourier_fit_w2_alt.order, fourier_fit_w2_alt.rmse, fourier_fit_w2_alt.r2);
    end
end
if processing_modes.s1
    fprintf('\nFourier parameters for W_1 (k, a_k, b_k, amplitude, phase):\n');
    disp(fourier_param_table_w1);
    if has_symmetric_peak_s1
        fprintf('Fourier parameters for W_1 at symmetric psi (k, a_k, b_k, amplitude, phase):\n');
        disp(fourier_param_table_w1_alt);
    end
end
if processing_modes.s2
    fprintf('Fourier parameters for W_2 (k, a_k, b_k, amplitude, phase):\n');
    disp(fourier_param_table_w2);
    if has_symmetric_peak_s2
        fprintf('Fourier parameters for W_2 at symmetric psi (k, a_k, b_k, amplitude, phase):\n');
        disp(fourier_param_table_w2_alt);
    end
end
if export_fourier_params
    if processing_modes.s1
        fprintf('Saved Fourier parameter CSV (W_1): %s\n', csv_path_w1);
    end
    if processing_modes.s2
        fprintf('Saved Fourier parameter CSV (W_2): %s\n', csv_path_w2);
    end
    fprintf('Saved Fourier parameter MAT: %s\n', mat_path_fit);
end
if processing_modes.s1
    fprintf('Gamma_1 at psi_plus^* (z=W_1) = %.12f\n', Gamma1_from_w(idx_max_s1));
    fprintf('Gamma_1 at psi_plus^* (z=cos) = %.12f\n', Gamma1_from_sine(idx_max_s1));
    if has_symmetric_peak_s1
        fprintf('Gamma_1 at symmetric psi_plus^* (z=W_1) = %.12f\n', Gamma1_from_w_alt(idx_max_s1_alt));
    end
end
if processing_modes.s2
    fprintf('Gamma_2 at psi_plus^* (z=W_2) = %.12f\n', Gamma2_from_w(idx_max_s2));
    fprintf('Gamma_2 at psi_plus^* (z=cos) = %.12f\n', Gamma2_from_sine(idx_max_s2));
    if has_symmetric_peak_s2
        fprintf('Gamma_2 at symmetric psi_plus^* (z=W_2) = %.12f\n', Gamma2_from_w_alt(idx_max_s2_alt));
    end
end

function [r_max, idx_max, idx_alt, has_symmetric_peak, r_alt] = find_symmetric_max_candidate(r_values, psi_grid)
    [r_max, idx_max] = max(r_values);
    idx_alt = NaN;
    has_symmetric_peak = false;
    r_alt = NaN;

    if numel(psi_grid) < 2 || ~isfinite(r_max)
        return;
    end

    psi_step = median(abs(diff(psi_grid)));
    if ~isfinite(psi_step) || psi_step <= 0
        psi_step = 1e-6;
    end

    tol_r = max(1e-12, 1e-9 * max(1, abs(r_max)));
    idx_max_all = find(abs(r_values - r_max) <= tol_r);
    if isempty(idx_max_all)
        return;
    end

    candidate_idx = idx_max_all(idx_max_all ~= idx_max);
    if isempty(candidate_idx)
        return;
    end

    psi_ref = psi_grid(idx_max);
    [~, nearest_idx] = min(abs(psi_grid(candidate_idx) + psi_ref));
    idx_alt = candidate_idx(nearest_idx);
    psi_alt = psi_grid(idx_alt);
    r_alt = r_values(idx_alt);

    has_symmetric_peak = abs(psi_ref + psi_alt) <= max(psi_step, 1e-8);
end

function w_values = compute_w_profile(phi_values, psi_plus, coeff, basis_types, m_order, n_order, orientation)
    switch orientation
        case 's1'
            w_minus = evaluate_exported_s(phi_values, phi_values - psi_plus, coeff, basis_types, m_order, n_order);
            w_plus = evaluate_exported_s(phi_values, phi_values + psi_plus, coeff, basis_types, m_order, n_order);
        case 's2'
            w_minus = evaluate_exported_s(phi_values - psi_plus, phi_values, coeff, basis_types, m_order, n_order);
            w_plus = evaluate_exported_s(phi_values + psi_plus, phi_values, coeff, basis_types, m_order, n_order);
        otherwise
            error('Unsupported orientation: %s', orientation);
    end
    w_values = w_minus - w_plus;
end

function gamma_values = compute_gamma_curve_from_z(psi_scan_values, phi_integral, z_values, coeff, basis_types, m_order, n_order, z_mean, orientation)
    gamma_values = nan(size(psi_scan_values));
    for i = 1:numel(psi_scan_values)
        psi_i = psi_scan_values(i);
        s_shifted = evaluate_shifted_surface(phi_integral, psi_i, coeff, basis_types, m_order, n_order, z_mean, orientation);
        gamma_values(i) = trapz(phi_integral, z_values .* s_shifted) / (2 * pi);
    end
end

function s_shifted = evaluate_shifted_surface(phi_values, psi_value, coeff, basis_types, m_order, n_order, z_mean, orientation)
    switch orientation
        case 's1'
            s_shifted = evaluate_exported_s(phi_values, phi_values - psi_value, coeff, basis_types, m_order, n_order) + z_mean;
        case 's2'
            s_shifted = evaluate_exported_s(phi_values - psi_value, phi_values, coeff, basis_types, m_order, n_order) + z_mean;
        otherwise
            error('Unsupported orientation: %s', orientation);
    end
end

function fit_data = get_fit_data_for_agent(gamma_export, agent_id, signal_role)
    agent_idx = find([gamma_export.agents.agent_id] == agent_id, 1, 'first');
    assert(~isempty(agent_idx), 'agent_id %d was not found in gamma_export.agents.', agent_id);

    switch lower(signal_role)
        case 'a2'
            fit_data = gamma_export.agents(agent_idx).a2_gamma;
        case 'derived'
            fit_data = gamma_export.agents(agent_idx).derived_gamma;
        otherwise
            error('signal_role must be ''a2'' or ''derived''.');
    end
end

function [coeff, basis_types, m_order, n_order, z_mean] = extract_basis_data(fit_data)
    assert(isfield(fit_data, 'coeff') && ~isempty(fit_data.coeff), 'No coeff was found in selected fit data.');
    assert(isfield(fit_data, 'basis_types') && ~isempty(fit_data.basis_types), 'No basis_types was found in selected fit data.');

    coeff = fit_data.coeff(:);
    basis_types = fit_data.basis_types;
    z_mean = 0;
    if isfield(fit_data, 'z_mean') && ~isempty(fit_data.z_mean) && isfinite(fit_data.z_mean(1))
        z_mean = double(fit_data.z_mean(1));
    else
        warning(['z_mean was not found in exported fit data. ', ...
            'Reconstructed surface baseline may differ from original plot. ', ...
            'Re-export gamma data with updated exporter to include z_mean.']);
    end

    m_order = [];
    n_order = [];
    if isfield(fit_data, 'basis_phi1_order')
        m_order = fit_data.basis_phi1_order(:);
    end
    if isfield(fit_data, 'basis_phi2_order')
        n_order = fit_data.basis_phi2_order(:);
    end

    if (isempty(m_order) || isempty(n_order)) && isfield(fit_data, 'basis_table') && istable(fit_data.basis_table)
        basis_table = fit_data.basis_table;
        if isempty(m_order) && ismember('phi1_order', basis_table.Properties.VariableNames)
            m_order = basis_table.phi1_order(:);
        end
        if isempty(n_order) && ismember('phi2_order', basis_table.Properties.VariableNames)
            n_order = basis_table.phi2_order(:);
        end
    end

    assert(~isempty(m_order) && ~isempty(n_order), 'No basis order arrays were found in selected fit data.');
    assert(numel(coeff) == numel(basis_types), 'Size mismatch: coeff and basis_types.');
    assert(numel(coeff) == numel(m_order), 'Size mismatch: coeff and basis_phi1_order.');
    assert(numel(coeff) == numel(n_order), 'Size mismatch: coeff and basis_phi2_order.');
end

function s_values = evaluate_exported_s(phi1, phi2, coeff, basis_types, m_order, n_order)
    s_values = zeros(size(phi1));

    for k = 1:numel(coeff)
        mk = m_order(k);
        nk = n_order(k);
        t = basis_types{k};

        switch t
            case 'constant'
                basis_value = ones(size(phi1));
            case 'phi1_cos'
                basis_value = cos(mk * phi1);
            case 'phi1_sin'
                basis_value = sin(mk * phi1);
            case 'phi2_cos'
                basis_value = cos(nk * phi2);
            case 'phi2_sin'
                basis_value = sin(nk * phi2);
            case 'mixed_cc'
                basis_value = cos(mk * phi1) .* cos(nk * phi2);
            case 'mixed_cs'
                basis_value = cos(mk * phi1) .* sin(nk * phi2);
            case 'mixed_sc'
                basis_value = sin(mk * phi1) .* cos(nk * phi2);
            case 'mixed_ss'
                basis_value = sin(mk * phi1) .* sin(nk * phi2);
            otherwise
                error('Unsupported basis type: %s', t);
        end

        s_values = s_values + coeff(k) * basis_value;
    end
end

function fit_result = fit_fourier_series_periodic(phi, y, max_order)
    phi_col = phi(:);
    y_col = y(:);
    max_order = max(0, round(double(max_order)));

    n_samples = numel(phi_col);
    design = zeros(n_samples, 1 + 2 * max_order);
    design(:, 1) = 1;

    for k = 1:max_order
        design(:, 2 * k) = cos(k * phi_col);
        design(:, 2 * k + 1) = sin(k * phi_col);
    end

    coeff_ls = design \ y_col;
    y_fit_col = design * coeff_ls;

    fit_result = struct();
    fit_result.order = max_order;
    fit_result.a0 = coeff_ls(1);
    fit_result.a = zeros(max_order, 1);
    fit_result.b = zeros(max_order, 1);
    if max_order > 0
        fit_result.a = coeff_ls(2:2:(2 * max_order));
        fit_result.b = coeff_ls(3:2:(2 * max_order + 1));
    end
    fit_result.y_fit = reshape(y_fit_col, size(phi));
    fit_result.rmse = sqrt(mean((y_col - y_fit_col) .^ 2));

    y_var = sum((y_col - mean(y_col)) .^ 2);
    fit_result.r2 = 1 - sum((y_col - y_fit_col) .^ 2) / max(y_var, eps);
end

function param_table = build_fourier_param_table(fit_result)
    k = (0:fit_result.order).';
    a_k = [fit_result.a0; fit_result.a(:)];
    b_k = [0; fit_result.b(:)];
    amplitude = sqrt(a_k .^ 2 + b_k .^ 2);
    phase = atan2(-b_k, a_k);

    param_table = table(k, a_k, b_k, amplitude, phase, ...
        'VariableNames', {'k', 'a_k', 'b_k', 'amplitude', 'phase'});
end
