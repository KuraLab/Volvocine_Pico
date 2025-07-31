% データ入力
alpha = [-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 -0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]*pi;
phase = [0 0 0.5*pi pi 0.9*pi 0.9*pi pi pi 0.9*pi 0.8*pi 0.85*pi -0.15*pi 0.1*pi 0 0 0.05*pi 0 0 0 0 0];
phase2 = [-0.129 -0.091 -0.032 -0.195 -0.099 -0.101 -0.054 -0.196 -0.307 -0.306 2.811 2.707 2.875 2.910 3.084 2.855 2.830 3.141 -1.490 -0.027 -0.129];

% プロット
figure;
plot(alpha, abs(phase2), '-o', 'LineWidth', 2);
xlabel('\alpha');
ylabel('\phi_2 - \phi_1');
xlim([-pi pi])
ylim([-pi pi])

% x軸をπベースで設定
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

% y軸をπベースで設定
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

grid on;
tuneFigure;
saveFigure;