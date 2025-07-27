% データ入力
alpha = [-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 -0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]*pi;
phase = [0 0 0.5*pi pi 0.9*pi 0.9*pi pi pi 0.9*pi 0.8*pi 0.85*pi -0.15*pi 0.1*pi 0 0 0.05*pi 0 0 0 0 0];

% プロット
figure;
plot(alpha/pi, phase/pi, '-o', 'LineWidth', 2);
xlabel('\alpha / \pi');
ylabel('位相差 / \pi');
grid on;
