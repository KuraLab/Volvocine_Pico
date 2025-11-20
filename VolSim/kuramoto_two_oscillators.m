% Simple Kuramoto model with two coupled oscillators and phase difference plot
% This script integrates the classic Kuramoto equations for two oscillators
% and visualizes how their phase difference evolves over time.

% Parameters
omega = pi*[2.50; 1.27];   % natural frequencies of the oscillators
K0 = 2;               % zeroth-harmonic coupling strength
K = -0.0;                % first-harmonic coupling strength
K2 =-0.5;              % second-harmonic coupling strength
K3 = 0.0;              % second-harmonic coupling strength

tspan = [0 200];         % simulation time span
tEval = linspace(tspan(1), tspan(2), 5000); % desired output times for smooth curves
numRuns = 10;
initialDiffs = (2*pi).*rand(1, numRuns) - pi; % random initial phase differences
theta1_initial = 0;     % reference phase for oscillator 1
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Integrate Kuramoto equations for each initial difference
oFun = @(t, theta) kuramoto_rhs(theta, omega, K0, K, K2, K3);
phaseDiffAll = zeros(numel(initialDiffs), numel(tEval));

for idx = 1:numel(initialDiffs)
    theta0 = [theta1_initial; wrap_to_pi(theta1_initial + initialDiffs(idx))];
    [t, theta] = ode45(oFun, tEval, theta0, options);
    phaseDiffAll(idx, :) = wrap_to_pi(2*theta(:, 2) - 1*theta(:, 1));
end

% Plot phase difference over time for all trajectories
figure;
hold on;
colors = lines(numel(initialDiffs));
for idx = 1:numel(initialDiffs)
    cleanSeries = break_phase_jumps(phaseDiffAll(idx, :));
    plot(t, cleanSeries, 'LineWidth', 1.2, 'Color', colors(idx, :));
end
hold off;
xlabel('time');
ylabel('$$\theta_1 - 2\theta_2$$','Interpreter','latex');
ylim([-pi pi]);
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
grid on;
if exist('tuneFigure', 'file')
    tuneFigure;
end

% Optional save
% saveFigure;

function dtheta = kuramoto_rhs(theta, omega, K0, K, K2, K3)
%KURAMOTO_RHS Right-hand side of the two-oscillator Kuramoto model.
%   theta: 2x1 vector of oscillator phases.
%   omega: 2x1 vector of natural frequencies.
%   K: scalar coupling strength for sin(diff).
%   K2: scalar coupling strength for sin(2*diff).

    phaseDiff = theta(2) - theta(1);
    coupling0 = (K0 / 2) * [cos(theta(1)); cos(theta(2))];
    coupling1 = (K / 2) * [sin(phaseDiff); sin(-phaseDiff)];
    coupling2 = (K2 / 2) * [sin(2*phaseDiff); sin(-2*phaseDiff)];
    coupling = coupling0 + coupling1 + coupling2;
    dtheta = omega + coupling;
end
function wrapped = wrap_to_pi(phi)
%WRAP_TO_PI Wrap angles to the interval [-pi, pi].

    wrapped = mod(phi + pi, 2*pi) - pi;
end

function cleanSeries = break_phase_jumps(series)
%BREAK_PHASE_JUMPS Insert NaNs when consecutive wrapped phases differ by >= pi.

    cleanSeries = series;
    jumpIdx = find(abs(diff(series)) >= pi);
    cleanSeries(jumpIdx + 1) = NaN;
end
