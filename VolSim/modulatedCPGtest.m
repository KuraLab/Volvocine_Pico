% forced_slave_sync.m
% Simulation: master oscillator with step-like (0/1) output forcing a slave oscillator
% Master: dtheta_m/dt = Omega_m
% Master output s(t) = 1 (on) or 0 (off) — step waveform depending on time/phase
% Slave (second-order default here): d^2 theta_s/dt^2 = (omega_s - dtheta_s/dt) + k * sin(theta_s + alpha) * s(t)
%   Here we add a linear term (omega_s - dtheta_s/dt) that pulls the velocity back to omega_s,
%   so that on average dtheta_s/dt ≈ omega_s and the system does not diverge.
%
% You can switch where the feedback is applied:
%   feedbackMode = 'acceleration'  -> feedback acts on d^2 theta_s/dt^2 (second-order model)
%   feedbackMode = 'velocity'      -> feedback acts directly on dtheta_s/dt (first-order model)
%
% Usage:
% Open this file in the MATLAB command window and simply run it.
% Example: >> forced_slave_sync

clear; close all;
% Parameters for master and slaves
Omega_m1 = 1.5*pi*0.5;   % Master angular frequency (rad/s) in first half
Omega_m2 = 1.5*pi*1.0;   % Master angular frequency (rad/s) in second half
numSlaves = 50;                              % Number of slaves
omega_s_vec = linspace(1.5*pi*0.4, 1.5*pi*1.1, numSlaves);  % Intrinsic slave frequencies
k_base = 5;            % Base feedback strength (mean over slaves)
k_vec = k_base * ones(size(omega_s_vec));  % Same k for all slaves
alpha = 0.0*pi;          % Phase offset
duty = 0.68;              % Master duty ratio (0~1). 0.5 means 50% duty

% Choose where the feedback term acts: 'acceleration' or 'velocity'
feedbackMode = getenv('FEEDBACK_MODE');
if isempty(feedbackMode)
    feedbackMode = 'velocity';  % default to velocity feedback
end

% --- Time discretization (Euler method) ---
tmax = 30;               % Simulation time (s)
dt   = 0.001;            % Time step (s)
t    = 0:dt:tmax;        % Uniform time grid
if size(t,1) == 1
    t = t.';             % Convert row vector to column vector if needed
end

% Moving-average settings for feedback term
fb_ma_T = 0.00001;                          % window length in time [s]
fb_ma_N = max(1, round(fb_ma_T/dt));        % window length in steps

% Moving-average settings for master input s(t)
s_ma_T = 0.5;                              % window length in time [s]
s_ma_N = max(1, round(s_ma_T/dt));          % window length in steps

% Initial conditions
theta_m0 = 0.0;          % Master initial phase
theta_s0 = linspace(0, 2*pi, numSlaves+1);  % Slave initial phases (distributed)
theta_s0 = theta_s0(1:end-1);
v_s0     = omega_s_vec;  % Slave initial angular velocity

N = numel(t);

% Pre-allocation
theta_m = zeros(N,1);   % Master phase
theta_s = zeros(N,numSlaves);   % Slave phases
v_s     = zeros(N,numSlaves);   % Slave angular velocities dtheta_s/dt
fb_raw  = zeros(N,numSlaves);   % raw feedback terms
fb_ma   = zeros(N,numSlaves);   % moving-averaged feedback terms
s_raw   = zeros(N,1);   % raw master output
svec    = zeros(N,1);   % filtered master output used for control and plotting

theta_m(1) = theta_m0;
theta_s(1,:) = theta_s0;
v_s(1,:)     = v_s0;

% Master period and step output (initialize with first-half frequency)
Omega_m = Omega_m1;
t_period = 2*pi / Omega_m;  % Physical period of the master (s)

% --- Time evolution via Euler method ---
for n = 1:N-1
    tn = t(n);

    % Switch master frequency at half of the simulation time
    if tn < tmax/2
        Omega_m = Omega_m1;
    else
        Omega_m = Omega_m2;
    end
    t_period = 2*pi / Omega_m;   % update master period accordingly

    % Master output s(t): 1 if sin(theta_m) >= -0.2, else 0
    s_t = double(sin(theta_m(n)) >= -0.3)-0.5;  % step output centered at 0 (0.5 or -0.5)
    s_raw(n) = s_t;
    i0s = max(1, n - s_ma_N + 1);
    s_t_in = mean(s_raw(i0s:n));
    svec(n) = s_t_in;   % store filtered signal used for control

    % Master: dtheta_m/dt = Omega_m
    dtheta_m = Omega_m;

    % --- Feedback term with moving average ---
    fb_pre = sin(theta_s(n,:) + alpha) * s_t_in;
    fb_raw(n,:) = k_vec .* fb_pre;
    i0 = max(1, n - fb_ma_N + 1);
    fb_ma(n,:) = mean(fb_raw(i0:n,:), 1);   % smoothed feedback term

    if strcmp(feedbackMode,'acceleration')
        % Slave: second-order equation
        %   v_s = dtheta_s/dt
        %   dv_s/dt = fb_ma(n)
        dv_s = fb_ma(n,:);

        % Forward Euler update (second-order)
        theta_m(n+1) = theta_m(n) + dt * dtheta_m;
        v_s(n+1,:)     = v_s(n,:)     + dt * dv_s;
        theta_s(n+1,:) = theta_s(n,:) + dt * v_s(n,:);  % Use v_s(n,:) (explicit Euler)

    else  % 'velocity' mode
        % Slave: first-order phase equation with smoothed feedback on velocity
        %   dtheta_s/dt = omega_s_vec + fb_ma(n,:)
        dtheta_s = omega_s_vec + fb_raw(n,:);

        % Forward Euler update (first-order)
        theta_m(n+1) = theta_m(n) + dt * dtheta_m;
        theta_s(n+1,:) = theta_s(n,:) + dt * dtheta_s;

        % Define an effective angular velocity for plotting
        v_s(n+1,:) = dtheta_s;   % instantaneous omega_s(t) in this mode
    end
end

% Set the last sample of svec equal to the last computed value
s_raw(end) = s_raw(end-1);
svec(end) = svec(end-1);

% Phase difference wrapped to [-pi, pi]
theta_m_mat = repmat(theta_m, 1, numSlaves);
phi = wrapToPi(theta_s - theta_m_mat);

% Collective Kuramoto order parameter among all slaves
z_slaves = exp(1i*theta_s);
z_order  = mean(z_slaves, 2);
R_slaves = abs(z_order);
phi_order_master = wrapToPi(angle(z_order) - theta_m);

% Prepare phase data for plotting without visual jumps (insert NaNs on big jumps)
phi_m = mod(theta_m, 2*pi);
phi_s = mod(theta_s, 2*pi);
z_order_phase = mod(angle(z_order), 2*pi);

thresh = pi;  % threshold for detecting jumps
jump_idx_m = find(abs(diff(phi_m)) >= thresh) + 1;
jump_idx_order = find(abs(diff(z_order_phase)) >= thresh) + 1;
jump_idx_s = cell(numSlaves,1);
for j = 1:numSlaves
    jump_idx_s{j} = find(abs(diff(phi_s(:,j))) >= thresh) + 1;
end

phi_m_plot = phi_m;
phi_s_plot = phi_s;
phi_order_plot = z_order_phase;
phi_m_plot(jump_idx_m) = NaN;
phi_order_plot(jump_idx_order) = NaN;
for j = 1:numSlaves
    phi_s_plot(jump_idx_s{j},j) = NaN;
end

% Visualization: master output and slave phases
figure('Name','Forced Multi-Slave Synchronization (Euler, selectable feedback)','NumberTitle','off','Position',[200 200 900 850]);

% Top: master output s(t)
subplot(2,1,1);
plot(t, svec, 'k','LineWidth',1.0);
ylim([-0.6 0.6]); ylabel('s(t)');

% Second: phases (with NaNs inserted to break jumps)
subplot(2,1,2);
h_slaves_phase = plot(t, phi_s_plot, 'r', 'LineWidth',0.5); hold on;
h_master_phase = plot(t, phi_m_plot, 'b', 'LineWidth',1.2);
h_order_phase = plot(t, phi_order_plot, 'g', 'LineWidth',1.2);
ylabel('phase (mod 2\pi)');
xlabel('time (s)');
ylim([0 2*pi]);
legend([h_slaves_phase(1), h_master_phase, h_order_phase], 'slaves', 'master', 'order phase'); 

% Simple synchronization check: variation of phase differences in last 10% of time
tf_final = t > (0.9 * tmax);
phi_std_final = std(phi(tf_final,:), 0, 1);
num_locked = sum(phi_std_final < 1e-2);
fprintf('Final phase-diff std (last 10%%): min %.4e / median %.4e / max %.4e rad\n', ...
    min(phi_std_final), median(phi_std_final), max(phi_std_final));
fprintf(' -> Locked slaves (std < 1e-2): %d / %d\n', num_locked, numSlaves);

% ===== Order parameters =====
% Collective Kuramoto order parameter among all slaves

% Per-slave order parameter relative to master (cumulative in time)
cum_complex = cumsum(exp(1i*phi), 1);
r_slave_t = abs(cum_complex ./ (1:N)');

figure('Name','Order Parameters (Multi-Slave)','NumberTitle','off','Position',[250 120 900 900]);

subplot(3,1,1);
plot(t, R_slaves, 'b', 'LineWidth',1.4);
ylim([0 1]); ylabel('R_{slaves}(t)');
title('Collective order parameter of all slaves');

subplot(3,1,2);
plot(t, r_slave_t, 'LineWidth',0.6);
ylim([0 1]); ylabel('r_j(t)');
title('Per-slave order parameter relative to master');

subplot(3,1,3);
plot(t, phi_order_master, 'm', 'LineWidth',1.2);
ylim([-pi pi]);
ylabel('\angle R_{slaves} - \theta_m');
xlabel('time (s)');
title('Phase difference: collective slave order phase vs master');

% ===== Complex-plane visualization on unit circle =====
% Master/slaves are projected as z = exp(i*theta), and the slave order
% parameter is shown both as the complex mean vector and its unit-phase point.
z_master = exp(1i * theta_m);                 % N x 1, on unit circle
z_order_u = exp(1i * angle(z_order));         % N x 1, phase-only on unit circle

figure('Name','Complex Unit-Circle View (Master/Slaves/Order)','NumberTitle','off', ...
    'Position',[1200 180 700 700]);
hFig = gcf;

th_c = linspace(0, 2*pi, 400);
plot(cos(th_c), sin(th_c), 'k--', 'LineWidth', 1.0); hold on;
axis equal; grid on;
xlim([-1.2 1.2]); ylim([-1.2 1.2]);
xlabel('Re'); ylabel('Im');
title('Complex-plane rotation on unit circle');

%h_master_tr = plot(NaN, NaN, 'b-', 'LineWidth', 1.0);
h_master    = plot(NaN, NaN, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
h_slaves    = plot(NaN, NaN, 'r.', 'MarkerSize', 15);
h_order_vec = plot([0 NaN], [0 NaN], 'g-', 'LineWidth', 5);
h_order     = plot(NaN, NaN, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
h_order_u   = plot(NaN, NaN, 'md', 'MarkerFaceColor', 'm', 'MarkerSize', 10);

legend({'unit circle','master','slaves','order vector','order (complex mean)','order phase on unit circle'}, ...
    'Location','southoutside');

% Animation/export settings
%saveVideo = ~strcmpi(getenv('SAVE_COMPLEX_VIDEO'), '0');  % 既定true, 環境変数で0なら無効
saveVideo = false;  % 既定true, 環境変数で0なら無効
videoFile = fullfile(pwd, 'complex_unit_circle.mp4');
videoFPS = 20;                         % output movie frame rate (lower is safer/lighter)
playbackSlowdown = 1.0;                % >1.0 slows on-screen playback (video speed unchanged)
maxVideoFrames = 300;                  % 保存フレーム上限（小さいほど高速）
captureDPI = 100;                      % printキャプチャ解像度（小さいほど高速）
maxDisplaySlaves = 40;                 % 描画するslave点数上限（計算自体は全slave）
captureMethod = getenv('VIDEO_CAPTURE_METHOD');
if isempty(captureMethod)
    captureMethod = 'print';           % 'print' (safer) or 'getframe'
end

anim_dt = 0.02;                                % visualization update interval [s]
anim_step_base = max(1, round(anim_dt / dt));
trail_T = 1.0;                                 % master trail length [s]
trail_N = max(2, round(trail_T / dt));

if numSlaves > maxDisplaySlaves
    slaveDrawIdx = round(linspace(1, numSlaves, maxDisplaySlaves));
else
    slaveDrawIdx = 1:numSlaves;
end

if saveVideo
    estFrames = ceil(N / anim_step_base);
    frameStride = max(1, ceil(estFrames / maxVideoFrames));
else
    frameStride = 1;
end
anim_step = anim_step_base * frameStride;
fprintf('Animation step: %d (base=%d, stride=%d), expected frames: %d\n', ...
    anim_step, anim_step_base, frameStride, ceil(N/anim_step));

if saveVideo
    vw = [];
    videoIsOpen = false;
    if strcmp(captureMethod, 'print')
        set(hFig, 'Renderer', 'painters');
    end
    try
        vw = VideoWriter(videoFile, 'MPEG-4');
        vw.FrameRate = videoFPS;
        open(vw);
        videoIsOpen = true;
        fprintf('Saving complex-plane animation (MPEG-4) to: %s\n', videoFile);
    catch
        [videoDir, videoName] = fileparts(videoFile);
        videoFile = fullfile(videoDir, [videoName '.avi']);
        try
            vw = VideoWriter(videoFile, 'Motion JPEG AVI');
            vw.FrameRate = videoFPS;
            open(vw);
            videoIsOpen = true;
            fprintf('MPEG-4 unavailable. Fallback to Motion JPEG AVI: %s\n', videoFile);
        catch ME
            saveVideo = false;
            fprintf('VideoWriter init failed. Continue without save: %s\n', ME.message);
        end
    end
end

for n = 1:anim_step:N
    i0 = max(1, n - trail_N + 1);

    %set(h_master_tr, 'XData', real(z_master(i0:n)), 'YData', imag(z_master(i0:n)));
    set(h_master,    'XData', real(z_master(n)),    'YData', imag(z_master(n)));
    set(h_slaves,    'XData', real(z_slaves(n,slaveDrawIdx)),  'YData', imag(z_slaves(n,slaveDrawIdx)));
    set(h_order_vec, 'XData', [0 real(z_order(n))], 'YData', [0 imag(z_order(n))]);
    set(h_order,     'XData', real(z_order(n)),     'YData', imag(z_order(n)));
    set(h_order_u,   'XData', real(z_order_u(n)),   'YData', imag(z_order_u(n)));

    title(sprintf('Complex-plane rotation on unit circle   t = %.2f s, |R| = %.3f', t(n), abs(z_order(n))));
    drawnow limitrate;

    if saveVideo && videoIsOpen
        try
            if strcmp(captureMethod, 'print')
                frameRGB = print(hFig, '-RGBImage', ['-r' num2str(captureDPI)]);
                writeVideo(vw, frameRGB);
            else
                writeVideo(vw, getframe(hFig));
            end
        catch ME
            fprintf('Frame write failed at n=%d. Stop video export: %s\n', n, ME.message);
            saveVideo = false;
            try
                close(vw);
            catch
            end
        end
    end

    if playbackSlowdown > 1.0
        pause(anim_dt * (playbackSlowdown - 1.0));
    end
end

if saveVideo && videoIsOpen
    close(vw);
    fprintf('Video saved: %s\n', videoFile);
elseif saveVideo && ~videoIsOpen
    fprintf('Video was not saved because VideoWriter could not be opened.\n');
end

function p = wrapToPi(x)
    % Wrap x into the range [-pi, pi]
    p = mod(x + pi, 2*pi) - pi;
end
