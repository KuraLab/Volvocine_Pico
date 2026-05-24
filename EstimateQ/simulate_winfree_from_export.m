function result = simulate_winfree_from_export(opts)
%SIMULATE_WINFREE_FROM_EXPORT Reconstruct s2 and z_opt from exported gamma data,
%then simulate the two-phase Winfree-type model:
%   dphi1/dt = omega1 + sigma * z(phi1) * s(phi1, phi2)
%   dphi2/dt = omega2 + sigma * z(phi2) * s(phi2, phi1)
%
% Usage:
%   result = simulate_winfree_from_export();
%   result = simulate_winfree_from_export(struct('omega', [240, 245], 'sigma', 0.02));
%
% Required export content (inside gamma_export):
%   - phase_agent_ids
%   - agents(:).a2_gamma or agents(:).derived_gamma
%   - coeff, basis_types, basis_phi1_order, basis_phi2_order (or basis_table)
%
% Output `result` fields:
%   - meta: simulation/export metadata
%   - t: time vector
%   - phi: [N x 2] phase trajectory
%   - s2_eval: function handle for reconstructed s2(phi1, phi2)
%   - z_opt_eval: function handle for reconstructed z_opt(theta)
%   - psi_plus_opt, scale_z_opt

if nargin < 1
    opts = struct();
end

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

% Defaults
opts = with_default(opts, 'mat_path', '');
opts = with_default(opts, 'signal_role', 'a2');   % 'a2' or 'derived'
opts = with_default(opts, 'agent_mode', 'phase_id_2'); % 'phase_id_2' | 'phase_id_1' | 'agent_id'
opts = with_default(opts, 'agent_id', NaN);
opts = with_default(opts, 'estimateq_root', script_dir);

opts = with_default(opts, 'omega', [240, 245]);
opts = with_default(opts, 'sigma', 0.02);
opts = with_default(opts, 'tspan', [0, 20]);
opts = with_default(opts, 'phi0', [0, pi/3]);

opts = with_default(opts, 'n_phi_integral', 2001);
opts = with_default(opts, 'n_psi_scan', 101);
opts = with_default(opts, 'psi_scan_min', -pi);
opts = with_default(opts, 'psi_scan_max', pi);
opts = with_default(opts, 'normalize_z_opt_power', true);

opts = with_default(opts, 'ode_rel_tol', 1e-8);
opts = with_default(opts, 'ode_abs_tol', 1e-10);

mat_path = resolve_mat_path(opts.mat_path, opts.estimateq_root);
S = load(mat_path, 'gamma_export');
assert(isfield(S, 'gamma_export'), 'The MAT file does not contain gamma_export.');
G = S.gamma_export;

assert(isfield(G, 'phase_agent_ids') && numel(G.phase_agent_ids) >= 2, ...
    'gamma_export.phase_agent_ids must contain at least two IDs.');
phase_ids = G.phase_agent_ids(:).';

switch lower(opts.agent_mode)
    case 'phase_id_2'
        agent_id = phase_ids(2);
    case 'phase_id_1'
        agent_id = phase_ids(1);
    case 'agent_id'
        assert(isfinite(opts.agent_id), 'opts.agent_id must be finite when agent_mode=''agent_id''.');
        agent_id = opts.agent_id;
    otherwise
        error('Unsupported opts.agent_mode: %s', opts.agent_mode);
end

fit_data = get_fit_data_for_agent(G, agent_id, opts.signal_role);
[coeff, basis_types, m_order, n_order, z_mean] = extract_basis_data(fit_data);

% Reconstruct z_opt from W2 profile at psi_plus maximizing R(psi_plus).
phi_integral = linspace(0, 2 * pi, opts.n_phi_integral);
psi_scan_values = linspace(opts.psi_scan_min, opts.psi_scan_max, opts.n_psi_scan);
R_values = nan(size(psi_scan_values));

for i = 1:numel(psi_scan_values)
    psi_i = psi_scan_values(i);
    w_i = compute_w_profile(phi_integral, psi_i, coeff, basis_types, m_order, n_order, 's2');
    R_values(i) = trapz(phi_integral, w_i .^ 2);
end

[~, idx_max] = max(R_values);
psi_plus_opt = psi_scan_values(idx_max);
z_opt_raw = compute_w_profile(phi_integral, psi_plus_opt, coeff, basis_types, m_order, n_order, 's2');

scale_z_opt = 1.0;
if opts.normalize_z_opt_power
    target_power = trapz(phi_integral, sin(phi_integral) .^ 2) / (2 * pi);
    z_power_raw = trapz(phi_integral, z_opt_raw .^ 2) / (2 * pi);
    scale_z_opt = sqrt(target_power / max(z_power_raw, eps));
end
z_opt_profile = scale_z_opt * z_opt_raw;

% Periodic interpolation for z_opt(theta)
z_opt_eval = @(theta) interp1_periodic(phi_integral, z_opt_profile, theta);

% Reconstructed s2(phi1, phi2) including baseline z_mean.
s2_eval = @(phi1, phi2) evaluate_exported_s(wrap_to_2pi(phi1), wrap_to_2pi(phi2), ...
    coeff, basis_types, m_order, n_order) + z_mean;

omega = opts.omega(:).';
assert(numel(omega) == 2, 'opts.omega must have two elements [omega1 omega2].');
sigma = opts.sigma;

odefun = @(t, y) [ ...
    omega(1) + sigma * z_opt_eval(y(1)) * s2_eval(y(1), y(2)); ...
    omega(2) + sigma * z_opt_eval(y(2)) * s2_eval(y(2), y(1)) ...
    ];

ode_opts = odeset('RelTol', opts.ode_rel_tol, 'AbsTol', opts.ode_abs_tol);
[t, phi] = ode45(odefun, opts.tspan, opts.phi0(:), ode_opts);

result = struct();
result.meta = struct();
result.meta.mat_path = mat_path;
result.meta.signal_role = opts.signal_role;
result.meta.agent_id = agent_id;
result.meta.omega = omega;
result.meta.sigma = sigma;
result.meta.tspan = opts.tspan;
result.meta.phi0 = opts.phi0;
result.meta.z_mean = z_mean;

result.t = t;
result.phi = phi;
result.psi_scan_values = psi_scan_values;
result.R_values = R_values;
result.psi_plus_opt = psi_plus_opt;
result.scale_z_opt = scale_z_opt;

result.z_opt_grid.theta = phi_integral;
result.z_opt_grid.value = z_opt_profile;

result.s2_eval = s2_eval;
result.z_opt_eval = z_opt_eval;
end

function mat_path = resolve_mat_path(mat_path, estimateq_root)
if ~isempty(mat_path) && ~isfile(mat_path)
    candidate_mat_path = fullfile(estimateq_root, mat_path);
    if isfile(candidate_mat_path)
        mat_path = candidate_mat_path;
    end
end

if isempty(mat_path)
    files = dir(fullfile(estimateq_root, 'Spring*', '*', 'gamma_exports', 'gamma_export_latest.mat'));
    assert(~isempty(files), 'No gamma_export_latest.mat was found under estimateq_root.');
    [~, newest_idx] = max([files.datenum]);
    mat_path = fullfile(files(newest_idx).folder, files(newest_idx).name);
end

assert(isfile(mat_path), 'MAT file was not found: %s', mat_path);
end

function v = with_default(s, field_name, default_value)
if ~isfield(s, field_name) || isempty(s.(field_name))
    s.(field_name) = default_value;
end
v = s;
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

function w_values = compute_w_profile(phi_values, psi_plus, coeff, basis_types, m_order, n_order, orientation)
switch orientation
    case 's2'
        w_minus = evaluate_exported_s(phi_values - psi_plus, phi_values, coeff, basis_types, m_order, n_order);
        w_plus = evaluate_exported_s(phi_values + psi_plus, phi_values, coeff, basis_types, m_order, n_order);
    case 's1'
        w_minus = evaluate_exported_s(phi_values, phi_values - psi_plus, coeff, basis_types, m_order, n_order);
        w_plus = evaluate_exported_s(phi_values, phi_values + psi_plus, coeff, basis_types, m_order, n_order);
    otherwise
        error('Unsupported orientation: %s', orientation);
end
w_values = w_minus - w_plus;
end

function yq = interp1_periodic(x, y, xq)
period = 2 * pi;
x0 = x(1);
xq_wrapped = mod(xq - x0, period) + x0;

% Remove duplicated endpoint when x spans a full period (e.g., linspace with both 0 and 2*pi).
x_col = x(:);
y_col = y(:);
if numel(x_col) >= 2 && abs((x_col(end) - x_col(1)) - period) <= 1e-12 * max(1, period)
    x_col = x_col(1:end-1);
    y_col = y_col(1:end-1);
end

x_ext = [x_col; x_col(1) + period];
y_ext = [y_col; y_col(1)];
yq = interp1(x_ext, y_ext, xq_wrapped, 'pchip');
end

function xw = wrap_to_2pi(x)
xw = mod(x, 2 * pi);
end
