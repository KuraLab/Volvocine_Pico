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
% Parameters for master and slave
Omega_m1 = 1.5*pi*1.0;   % Master angular frequency (rad/s) in first half
Omega_m2 = 1.5*pi*0.7;   % Master angular frequency (rad/s) in second half
omega_s  = 1.5*pi*1.0;   % Target slave angular frequency (rad/s)
k = 10.0;                 % Feedback strength
alpha = 0.0*pi;          % Phase offset
duty = 0.68;              % Master duty ratio (0~1). 0.5 means 50% duty

% Choose where the feedback term acts: 'acceleration' or 'velocity'
feedbackMode = 'acceleration';   % <--- change here if you want acceleration feedback
%feedbackMode = 'velocity';        % <--- change here if you want velocity feedback

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

% Initial conditions
theta_m0 = 0.0;          % Master initial phase
theta_s0 = 0.5;          % Slave initial phase
v_s0     = omega_s;      % Slave initial angular velocity (start from omega_s)

N = numel(t);

% Pre-allocation
theta_m = zeros(N,1);   % Master phase
theta_s = zeros(N,1);   % Slave phase
v_s     = zeros(N,1);   % Slave angular velocity dtheta_s/dt
fb_raw  = zeros(N,1);   % raw feedback term k*sin(theta_s+alpha)*s_t
fb_ma   = zeros(N,1);   % moving-averaged feedback term
svec    = zeros(N,1);   % master output used for control and plotting

theta_m(1) = theta_m0;
theta_s(1) = theta_s0;
v_s(1)     = v_s0;

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

    % Master output s(t): time-based step wave (centered around 0)
    s_t = double( mod(tn, t_period) < duty * t_period ) - 0.5;
    svec(n) = s_t;   % store the exact signal used for control

    % Master: dtheta_m/dt = Omega_m
    dtheta_m = Omega_m;

    % --- Feedback term with moving average ---
    fb_raw(n) = k * sin(theta_s(n) + alpha) * s_t;
    i0 = max(1, n - fb_ma_N + 1);
    fb_ma(n) = mean(fb_raw(i0:n));   % smoothed feedback term

    if strcmp(feedbackMode,'acceleration')
        % Slave: second-order equation
        %   v_s = dtheta_s/dt
        %   dv_s/dt = (omega_s - v_s) + fb_ma(n)
        %dv_s = (omega_s - v_s(n)) + fb_ma(n);
        dv_s = fb_ma(n);

        % Forward Euler update (second-order)
        theta_m(n+1) = theta_m(n) + dt * dtheta_m;
        v_s(n+1)     = v_s(n)     + dt * dv_s;
        theta_s(n+1) = theta_s(n) + dt * v_s(n);  % Use v_s(n) (explicit Euler)

    else  % 'velocity' mode
        % Slave: first-order phase equation with smoothed feedback on velocity
        %   dtheta_s/dt = omega_s + fb_ma(n)
        dtheta_s = omega_s + fb_ma(n);

        % Forward Euler update (first-order)
        theta_m(n+1) = theta_m(n) + dt * dtheta_m;
        theta_s(n+1) = theta_s(n) + dt * dtheta_s;

        % Define an effective angular velocity for plotting
        v_s(n+1) = dtheta_s;   % instantaneous omega_s(t) in this mode
    end
end

% Set the last sample of svec equal to the last computed value
svec(end) = svec(end-1);

% Phase difference wrapped to [-pi, pi]
phi = wrapToPi(theta_s - theta_m);

% Prepare phase data for plotting without visual jumps (insert NaNs on big jumps)
phi_m = mod(theta_m, 2*pi);
phi_s = mod(theta_s, 2*pi);

thresh = pi;  % threshold for detecting jumps
jump_idx_m = find(abs(diff(phi_m)) >= thresh) + 1;
jump_idx_s = find(abs(diff(phi_s)) >= thresh) + 1;

phi_m_plot = phi_m;
phi_s_plot = phi_s;
phi_m_plot(jump_idx_m) = NaN;
phi_s_plot(jump_idx_s) = NaN;

% Also break the phase-difference plot on large jumps
phi_plot = phi;
jump_idx_phi = find(abs(diff(phi)) >= thresh) + 1;
phi_plot(jump_idx_phi) = NaN;

% Visualization: master output, slave phase, phase difference, and angular velocity
figure('Name','Forced Slave Synchronization (Euler, selectable feedback)','NumberTitle','off','Position',[200 200 900 800]);

% Top: master output s(t)
subplot(4,1,1);
plot(t, svec, 'k','LineWidth',1.0);
ylim([-0.6 0.6]); ylabel('s(t)');

% Second: phases (with NaNs inserted to break jumps)
subplot(4,1,2);
plot(t, phi_m_plot, 'b', 'LineWidth',1.2); hold on;
plot(t, phi_s_plot, 'r', 'LineWidth',1.0);
ylabel('phase (mod 2\pi)'); legend('master','slave'); 

% Third: phase difference (also broken on jumps)
subplot(4,1,3);
plot(t, phi_plot, 'm','LineWidth',1.2);
ylim([-pi pi]); ylabel('\phi = \theta_s - \theta_m');
% Bottom: slave angular velocity
subplot(4,1,4);
plot(t, v_s, 'k','LineWidth',1.0);
ylim([0, max(v_s)]);              % fix lower bound at 0
ylabel('\omega_s (rad/s)'); xlabel('time (s)');

% Simple synchronization check: small variation of phase difference in last 10% of time
tf_final = t > (0.9 * tmax);
phi_std_final = std(phi(tf_final));
fprintf('Final phase-diff std (last 10%%): %.4e rad\n', phi_std_final);
if phi_std_final < 1e-2
    fprintf(' -> Converged: likely phase-locked (small variation of phase difference)\n');
else
    fprintf(' -> Not converged: phase difference varies (possibly not synchronized)\n');
end

% ----- Local functions ----

function s = master_output(theta_m_vals, duty_cycle)
    % theta_m_vals can be scalar or vector
    % Output is a 0/1 array of the same size
    ph = mod(theta_m_vals, 2*pi);
    s = double(ph < (2*pi*duty_cycle));
end

function p = wrapToPi(x)
    % Wrap x into the range [-pi, pi]
    p = mod(x + pi, 2*pi) - pi;
end

% ===== Optional: example for parameter sweep (commented out) =====
%{
% If you want to perform a parameter sweep with the Euler version,
% it is convenient to wrap the above loop into a function and call it from here.
%}
