% forced_slave_sync.m
% シミュレーション: マスターがステップ型(0/1)出力を出し、スレーブに強制を与える系
% マスター: dtheta_m/dt = Omega_m
% マスター出力 s(t) = 1 (on) または 0 (off) — step型、位相に依存
% スレーブ: dtheta_s/dt = omega_s + k * sin(theta_s + alpha) * s(t)
%
% 使い方:
% MATLAB のコマンドウィンドウでこのファイルを開き、実行するだけです。
% 例: >> forced_slave_sync

clear; close all;
% パラメータ
Omega_m = 2*pi*1.0;    % マスターの角周波数 (rad/s) — 1 Hz
omega_s  = 2*pi*0.95;   % スレーブの固有角周波数 (rad/s)
k = 1.0;                % フィードバック強度
alpha = 0.0;            % 位相オフセット
duty = 0.5;             % マスターのオン比 (0~1)。0.5 は 50% duty

tmax = 100;             % シミュレーション時間 (s)
tspan = [0 tmax];

% 初期位相 (rad)
theta_m0 = 0.0;
theta_s0 = 0.5;
y0 = [theta_m0; theta_s0];

opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

% ODE を解く
[t,y] = ode45(@odefun, tspan, y0, opts);

theta_m = y(:,1);
theta_s = y(:,2);

% マスター出力 s(t) を計算
svec = master_output(theta_m, duty);

% 位相差を -pi..pi にラップ
phi = wrapToPi(theta_s - theta_m);

% 可視化
figure('Name','Forced Slave Synchronization','NumberTitle','off','Position',[200 200 900 700]);

subplot(3,1,1);
plot(t, svec, 'k','LineWidth',1.2);
ylim([-0.1 1.1]); ylabel('s(t)'); title('Master output (step 0/1)');

subplot(3,1,2);
plot(t, mod(theta_m,2*pi), 'b', 'LineWidth',1.2); hold on;
plot(t, mod(theta_s,2*pi), 'r', 'LineWidth',1.0);
ylabel('phase (mod 2\pi)'); legend('master','slave'); title('Phases (mod 2\pi)');

subplot(3,1,3);
plot(t, phi, 'm','LineWidth',1.2);
ylim([-pi pi]); ylabel('\phi = \theta_s - \theta_m'); xlabel('time (s)');
title('Phase difference (wrapped to [-\pi,\pi])');

% 簡単な同期判定: 最後の 10% 時間での位相差の変動が小さいか
tf_final = t > (0.9 * tmax);
phi_std_final = std(phi(tf_final));
fprintf('Final phase-diff std (last 10%%): %.4e rad\n', phi_std_final);
if phi_std_final < 1e-2
    fprintf(' -> 収束: 位相ロックの可能性が高い (位相差の変動小)
');
else
    fprintf(' -> 非収束: 位相差が変動 (同期していない可能性)
');
end

% ----- ローカル関数 -----
function dy = odefun(~, y)
    theta_m = y(1);
    theta_s = y(2);
    s = master_output(theta_m, duty);
    dtheta_m = Omega_m;
    dtheta_s = omega_s + k * sin(theta_s + alpha) * s;
    dy = [dtheta_m; dtheta_s];
end

function s = master_output(theta_m_vals, duty_cycle)
    % theta_m_vals はスカラーまたはベクトル
    % 出力は同じサイズの 0/1 配列
    ph = mod(theta_m_vals, 2*pi);
    s = double(ph < (2*pi*duty_cycle));
end

function p = wrapToPi(x)
    % x を -pi..pi の範囲にラップ
    p = mod(x + pi, 2*pi) - pi;
end

% ===== オプション: パラメータスイープの例（コメント解除して実行） =====
%{ 
% k を増やして同期のしきい値を探す簡単なスイープ
% ks = linspace(0,3,15);
% final_phi_std = zeros(size(ks));
% for ii=1:length(ks)
%     k_tmp = ks(ii);
%     % 一時的な odefun を用意
%     odefun_tmp = @(t,y) [Omega_m; omega_s + k_tmp*sin(y(2)+alpha)*master_output(y(1),duty)];
%     [tt, yy] = ode45(odefun_tmp, tspan, y0, opts);
%     phi_tmp = wrapToPi(yy(:,2) - yy(:,1));
%     final_phi_std(ii) = std(phi_tmp(tt > (0.9*tmax)));
% end
% figure; plot(ks, final_phi_std, '-o'); xlabel('k'); ylabel('final phase std');
% title('Parameter sweep: dependence on k');
%}
