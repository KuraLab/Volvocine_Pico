% WinfreeArnold
% Masterの固有角速度を一定に固定し、k と omega_s を掃引して
% 平均周波数比 <omega_slave>/<omega_master> のカラーマップを作る。

close all; clc;
fprintf('Running script: %s\n', mfilename('fullpath'));

% ===== 固定パラメータ =====
Omega_m = 1.5*pi*0.5;    % Master angular frequency (rad/s) [固定]
alpha = 0.0*pi;          % Phase offset
k2 = 1.0;                % 2nd harmonic feedback gain
beta = 0.0*pi;           % 2nd harmonic phase offset
duty = 0.68;             % Master duty ratio (0~1)

% beta 掃引（スライダー表示用）
enableBetaSlider = true;
beta_min = -pi;
beta_max = pi;
beta_count = 41;

% ===== 実行プリセット =====
% 'fast' : まず形を見る用（高速）
% 'full' : 高解像度（低速）
runMode = 'fast';

% ===== 掃引範囲の手動指定 =====
% true にすると下の min/max/step 指定を使用（runModeより優先）
useCustomSweep = true;

% k の範囲
k_min = 0;
k_max = 5;
k_step = 0.1;

% omega_s の範囲（Omega_m比で指定）
omega_ratio_min = 0.4;   % omega_s / Omega_m
omega_ratio_max = 2.4;   % omega_s / Omega_m
omega_ratio_step = 0.05; % 刻み

if useCustomSweep
    if k_step <= 0 || omega_ratio_step <= 0
        error('k_step and omega_ratio_step must be positive.');
    end
    if k_max < k_min || omega_ratio_max < omega_ratio_min
        error('max must be >= min for both k and omega ratio.');
    end

    k_vals = k_min:k_step:k_max;
    omega_s_vals = (omega_ratio_min:omega_ratio_step:omega_ratio_max) * Omega_m;

    % 時間設定は runMode で選ぶ
    switch runMode
        case 'fast'
            tmax = 15;
            dt = 0.002;
        case 'full'
            tmax = 30;
            dt = 0.001;
        otherwise
            error('Unknown runMode: %s', runMode);
    end
else
    switch runMode
        case 'fast'
            k_vals = linspace(0, 20, 41);
            omega_s_vals = linspace(0.7*Omega_m, 1.3*Omega_m, 61);
            tmax = 15;
            dt = 0.002;
        case 'full'
            k_vals = linspace(0, 20, 81);
            omega_s_vals = linspace(0.6*Omega_m, 1.4*Omega_m, 121);
            tmax = 30;
            dt = 0.001;
        otherwise
            error('Unknown runMode: %s', runMode);
    end
end

% ===== 時間離散 =====
t = (0:dt:tmax).';

% 初期条件
theta_s0 = 0.5;

% 解析区間（過渡を捨てる）
analysis_ratio = 0.5;  % 後半50%で平均周波数を評価

N = numel(t);
n0 = max(1, floor((1 - analysis_ratio) * N));
t_period = 2*pi / Omega_m;
svec = double(mod(t(1:end-1), t_period) < duty * t_period) - 0.5;
dtheta_m_window = Omega_m * (t(end) - t(n0));

if enableBetaSlider
    beta_vals = linspace(beta_min, beta_max, beta_count);
else
    beta_vals = beta;
end

% 結果配列: 行=k, 列=omega_s, 3次元目=beta index
freq_ratio_cube = zeros(numel(k_vals), numel(omega_s_vals), numel(beta_vals));

% ===== 掃引計算 =====
tic;
for ib = 1:numel(beta_vals)
    beta_now = beta_vals(ib);
    for ik = 1:numel(k_vals)
        k = k_vals(ik);
        for io = 1:numel(omega_s_vals)
            omega_s = omega_s_vals(io);
            freq_ratio_cube(ik, io, ib) = simulate_avg_freq_ratio( ...
                omega_s, k, k2, alpha, beta_now, dt, ...
                theta_s0, N, n0, svec, dtheta_m_window);
        end
    end
    fprintf('Beta progress: %d/%d (beta=%.3f rad), elapsed %.1f s\n', ...
        ib, numel(beta_vals), beta_now, toc);
end

freq_ratio_map = freq_ratio_cube(:, :, 1);

% ===== 可視化（アーノルドの舌） =====
figure('Name','Arnold Tongue: <omega_slave>/<omega_master>', ...
       'NumberTitle','off', 'Position',[200 180 920 640]);

hImg = imagesc(omega_s_vals / Omega_m, k_vals, freq_ratio_map);
axis xy;
xlabel('\omega_s / \Omega_m');
ylabel('k');
tTitle = title(sprintf('Average frequency ratio map: <\\omega_{slave}> / <\\omega_{master}> (\\beta = %.3f rad)', beta_vals(1)));
cb = colorbar;
ylabel(cb, '<\omega_{slave}> / <\omega_{master}>');
colormap(turbo);
caxis([min(freq_ratio_cube(:)), max(freq_ratio_cube(:))]);

if enableBetaSlider && numel(beta_vals) > 1
    fig = gcf;
    tBeta = uicontrol(fig, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.15, 0.02, 0.7, 0.035], ...
        'String', sprintf('beta = %.3f rad (%.1f deg)', beta_vals(1), rad2deg(beta_vals(1))), ...
        'HorizontalAlignment', 'center');

    sBeta = uicontrol(fig, 'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.15, 0.06, 0.7, 0.03], ...
        'Min', 1, 'Max', numel(beta_vals), 'Value', 1, ...
        'SliderStep', [1/(numel(beta_vals)-1), 5/(numel(beta_vals)-1)]);

    setappdata(fig, 'hImg', hImg);
    setappdata(fig, 'tTitle', tTitle);
    setappdata(fig, 'tBeta', tBeta);
    setappdata(fig, 'freq_ratio_cube', freq_ratio_cube);
    setappdata(fig, 'beta_vals', beta_vals);
    set(sBeta, 'Callback', @on_beta_slider_changed);
end

fprintf('Sweep done: %d x %d x %d points\n', numel(k_vals), numel(omega_s_vals), numel(beta_vals));

% ----- Local function -----

function ratio = simulate_avg_freq_ratio(omega_s, k, k2, alpha, beta, dt, ...
                                         theta_s0, N, n0, svec, dtheta_m_window)
    theta_s = theta_s0;
    theta_s_n0 = theta_s0;

    for n = 1:N-1
        if n == n0
            theta_s_n0 = theta_s;
        end

        fb_term = k * (sin(theta_s + alpha) + sin(0.5*theta_s + beta)) * svec(n);

        dtheta_s = omega_s + fb_term;
        theta_s = theta_s + dt * dtheta_s;
    end

    dtheta_s_window = theta_s - theta_s_n0;
    ratio = dtheta_s_window / dtheta_m_window;
end

function on_beta_slider_changed(src, ~)
    fig = ancestor(src, 'figure');

    idx = round(get(src, 'Value'));
    set(src, 'Value', idx);

    hImg = getappdata(fig, 'hImg');
    tTitle = getappdata(fig, 'tTitle');
    tBeta = getappdata(fig, 'tBeta');
    freq_ratio_cube = getappdata(fig, 'freq_ratio_cube');
    beta_vals = getappdata(fig, 'beta_vals');

    beta_now = beta_vals(idx);
    set(hImg, 'CData', freq_ratio_cube(:, :, idx));
    set(tTitle, 'String', sprintf('Average frequency ratio map: <\\omega_{slave}> / <\\omega_{master}> (\\beta = %.3f rad)', beta_now));
    set(tBeta, 'String', sprintf('beta = %.3f rad (%.1f deg)', beta_now, rad2deg(beta_now)));
end
