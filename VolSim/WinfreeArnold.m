% WinfreeArnold
% Masterの固有角速度を一定に固定し、k と omega_s を掃引して
% 平均周波数比 <omega_slave>/<omega_master> のカラーマップを作る。
% 矩形入力（カクカク）と正弦入力（sin）を同条件で比較する。

close all; clc;
fprintf('Running script: %s\n', mfilename('fullpath'));

% ===== 固定パラメータ =====
Omega_m = 1.5*pi*0.5;    % Master angular frequency (rad/s) [固定]
alpha = 0.0*pi;          % Phase offset
k2 = 0.0;                % 2nd harmonic feedback gain
beta = -0.25*pi;           % 2nd harmonic phase offset
duty_50 = 0.5;           % Master duty ratio (0~1)
duty_70 = 0.7;           % Master duty ratio (0~1)

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
k_step = 0.2;

% omega_s の範囲（Omega_m比で指定）
omega_ratio_min = 0.0;   % omega_s / Omega_m
omega_ratio_max = 2.4;   % omega_s / Omega_m
omega_ratio_step = 0.02; % 刻み

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
svec_square_50 = double(mod(t(1:end-1), t_period) < duty_50 * t_period) - 0.5;
svec_square_70 = double(mod(t(1:end-1), t_period) < duty_70 * t_period) - 0.5;
svec_sine = 0.5 * sin(Omega_m * t(1:end-1));
dtheta_m_window = Omega_m * (t(end) - t(n0));

% 結果行列: 行=k, 列=omega_s
freq_ratio_map_square_50 = zeros(numel(k_vals), numel(omega_s_vals));
freq_ratio_map_square_70 = zeros(numel(k_vals), numel(omega_s_vals));
freq_ratio_map_sine = zeros(numel(k_vals), numel(omega_s_vals));

% ===== 掃引計算 =====
tic;
for ik = 1:numel(k_vals)
    k = k_vals(ik);
    for io = 1:numel(omega_s_vals)
        omega_s = omega_s_vals(io);
        freq_ratio_map_square_50(ik, io) = simulate_avg_freq_ratio( ...
            omega_s, k, k2, alpha, beta, dt, ...
            theta_s0, N, n0, svec_square_50, dtheta_m_window);

        freq_ratio_map_square_70(ik, io) = simulate_avg_freq_ratio( ...
            omega_s, k, k2, alpha, beta, dt, ...
            theta_s0, N, n0, svec_square_70, dtheta_m_window);

        freq_ratio_map_sine(ik, io) = simulate_avg_freq_ratio( ...
            omega_s, k, k2, alpha, beta, dt, ...
            theta_s0, N, n0, svec_sine, dtheta_m_window);
    end
    if mod(ik, max(1, round(numel(k_vals)/10))) == 0
        fprintf('Progress: %d/%d rows (%.0f%%), elapsed %.1f s\n', ...
            ik, numel(k_vals), 100*ik/numel(k_vals), toc);
    end
end

% ===== 入力波形プロット =====
show_periods = 3;
t_show_max = show_periods * t_period;
idx_show = t(1:end-1) <= t_show_max;

% Square duty=0.5
figure('Name','Master input: square duty=0.5', ...
       'NumberTitle','off', 'Position',[80 480 600 340]);
plot(t(idx_show), svec_square_50(idx_show), 'LineWidth', 1.6);
grid on; ylim([-0.7, 0.7]); xlabel('t [s]'); ylabel('s(t)');
tuneFigure;
%saveFigure;

% Square duty=0.7
figure('Name','Master input: square duty=0.7', ...
       'NumberTitle','off', 'Position',[700 480 600 340]);
plot(t(idx_show), svec_square_70(idx_show), 'LineWidth', 1.6);
grid on; ylim([-0.7, 0.7]); xlabel('t [s]'); ylabel('s(t)');
tuneFigure;
%saveFigure;

% Sine
figure('Name','Master input: sine', ...
       'NumberTitle','off', 'Position',[1320 480 600 340]);
plot(t(idx_show), svec_sine(idx_show), 'LineWidth', 1.6);
grid on; ylim([-0.7, 0.7]); xlabel('t [s]'); ylabel('s(t)');
tuneFigure;
%saveFigure;

% ===== フーリエ級数展開（各高調波の寄与） =====
n_harmonics = 10;
% 1周期分のデータで計算
idx_1period = t(1:end-1) <= t_period;
t_1period = t(idx_1period);

% Square duty=0.5
fourier_amps_50 = compute_fourier_amplitudes(t_1period, svec_square_50(idx_1period), Omega_m, n_harmonics);
figure('Name','Fourier harmonics: square duty=0.5', ...
       'NumberTitle','off', 'Position',[80 890 600 340]);
bar(0:n_harmonics, fourier_amps_50, 'FaceColor', [0.2 0.4 0.8]);
grid on;
xlabel('Harmonic number'); ylabel('Amplitude');
xlim([-0.5, n_harmonics+0.5]);
tuneFigure;
%saveFigure;

% Square duty=0.7
fourier_amps_70 = compute_fourier_amplitudes(t_1period, svec_square_70(idx_1period), Omega_m, n_harmonics);
figure('Name','Fourier harmonics: square duty=0.7', ...
       'NumberTitle','off', 'Position',[700 890 600 340]);
bar(0:n_harmonics, fourier_amps_70, 'FaceColor', [0.8 0.4 0.2]);
grid on;
xlabel('Harmonic number'); ylabel('Amplitude');
xlim([-0.5, n_harmonics+0.5]);
tuneFigure;
%saveFigure; 

% Sine
fourier_amps_sine = compute_fourier_amplitudes(t_1period, svec_sine(idx_1period), Omega_m, n_harmonics);
figure('Name','Fourier harmonics: sine', ...
       'NumberTitle','off', 'Position',[1320 890 600 340]);
bar(0:n_harmonics, fourier_amps_sine, 'FaceColor', [0.2 0.7 0.4]);
grid on;
xlabel('Harmonic number'); ylabel('Amplitude');
xlim([-0.5, n_harmonics+0.5]);
tuneFigure;
%saveFigure;

% ===== アーノルド舌プロット =====
% Square duty=0.5
figure('Name','Arnold Tongue: square duty=0.5', ...
    'NumberTitle','off', 'Position',[80 50 600 400]);
imagesc(omega_s_vals / Omega_m, k_vals, freq_ratio_map_square_50);
axis xy;
xlabel('$$\omega / \Omega$$'); ylabel('k');
cb1 = colorbar;
ylabel(cb1, '$$\bar{\dot{\phi}} / \Omega$$', 'Interpreter', 'latex');
colormap(gca, turbo);
tuneFigure;
%saveFigure;

% Square duty=0.7
figure('Name','Arnold Tongue: square duty=0.7', ...
    'NumberTitle','off', 'Position',[700 50 600 400]);
imagesc(omega_s_vals / Omega_m, k_vals, freq_ratio_map_square_70);
axis xy;
xlabel('$$\omega / \Omega$$'); ylabel('k');
cb2 = colorbar;
ylabel(cb2, '$$\bar{\dot{\phi}} / \Omega$$', 'Interpreter', 'latex');
colormap(gca, turbo);
tuneFigure;
%saveFigure;

% Sine
figure('Name','Arnold Tongue: sine input', ...
    'NumberTitle','off', 'Position',[1320 50 600 400]);
imagesc(omega_s_vals / Omega_m, k_vals, freq_ratio_map_sine);
axis xy;
xlabel('$$\omega / \Omega$$'); ylabel('k');
cb3 = colorbar;
ylabel(cb3, '$$\bar{\dot{\phi}} / \Omega$$', 'Interpreter', 'latex');
colormap(gca, turbo);
tuneFigure;
%saveFigure;

fprintf('Sweep done: %d x %d points (3 conditions)\n', numel(k_vals), numel(omega_s_vals));

% ----- Local function -----

function ratio = simulate_avg_freq_ratio(omega_s, k, k2, alpha, beta, dt, ...
                                         theta_s0, N, n0, svec, dtheta_m_window)
    theta_s = theta_s0;
    theta_s_n0 = theta_s0;

    for n = 1:N-1
        if n == n0
            theta_s_n0 = theta_s;
        end

        fb_term = k * (sin(theta_s + alpha) + k2 * sin(2*theta_s + beta)) * svec(n);

        dtheta_s = omega_s + fb_term;
        theta_s = theta_s + dt * dtheta_s;
    end

    dtheta_s_window = theta_s - theta_s_n0;
    ratio = dtheta_s_window / dtheta_m_window;
end

function amps = compute_fourier_amplitudes(t, s, omega0, n_max)
    % t: 時刻ベクトル (1周期分)
    % s: 信号ベクトル
    % omega0: 基本角周波数
    % n_max: 計算する最大次数
    
    T = 2*pi / omega0;  % 周期
    amps = zeros(n_max+1, 1);
    
    % DC成分 (n=0)
    a0 = trapz(t, s) / T;
    amps(1) = abs(a0);
    
    % 高調波成分 (n=1,2,...)
    for n = 1:n_max
        an = (2/T) * trapz(t, s .* cos(n * omega0 * t));
        bn = (2/T) * trapz(t, s .* sin(n * omega0 * t));
        amps(n+1) = sqrt(an^2 + bn^2);
    end
end
