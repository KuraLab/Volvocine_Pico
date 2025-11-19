% Simple Kuramoto model with two coupled oscillators and phase difference plot
% This script integrates the classic Kuramoto equations for two oscillators
% and visualizes how their phase difference evolves over time.

% Parameters
omega = [1.00; 1.00];   % natural frequencies of the oscillators
K = 0.1;                % first-harmonic coupling strength
K2 = 0.05;              % second-harmonic coupling strength

tspan = [0 80];         % simulation time span
tEval = linspace(tspan(1), tspan(2), 5000); % desired output times for smooth curves
theta0 = (2*pi).*rand(2,1) - pi;  % random initial phases in [-pi, pi]
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Integrate Kuramoto equations
oFun = @(t, theta) kuramoto_rhs(theta, omega, K, K2);
[t, theta] = ode45(oFun, tEval, theta0, options);

% Phase difference wrapped to [-pi, pi]
phaseDiff = wrap_to_pi(theta(:, 1) - theta(:, 2));

% Plot phase difference over time
figure;
plot(t, phaseDiff, 'LineWidth', 1.5);
xlabel('time');
ylabel('$$\Delta\theta = \theta_1 - \theta_2$$','Interpreter','latex');
ylim([-pi pi]);
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
grid on;
if exist('tuneFigure', 'file')
    tuneFigure;
end

% Optional save
% saveFigure;

function dtheta = kuramoto_rhs(theta, omega, K, K2)
%KURAMOTO_RHS Right-hand side of the two-oscillator Kuramoto model.
%   theta: 2x1 vector of oscillator phases.
%   omega: 2x1 vector of natural frequencies.
%   K: scalar coupling strength for sin(diff).
%   K2: scalar coupling strength for sin(2*diff).

    phaseDiff = theta(2) - theta(1);
    coupling1 = (K / 2) * [sin(phaseDiff); -sin(phaseDiff)];
    coupling2 = (K2 / 2) * [sin(2*phaseDiff); -sin(2*phaseDiff)];
    coupling = coupling1 + coupling2;
    dtheta = omega + coupling;
end

function wrapped = wrap_to_pi(phi)
%WRAP_TO_PI Wrap angles to the interval [-pi, pi].

    wrapped = mod(phi + pi, 2*pi) - pi;
end
