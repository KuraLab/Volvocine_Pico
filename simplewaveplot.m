% phi の範囲を -pi から pi まで取得
phi = linspace(-pi, pi, 1000);

% 関数値
y = 0.3*sin(phi) + 0.9*sin(2*phi);

% プロット
figure;
plot(phi, y, 'LineWidth', 1.5);
xlabel('$$\psi$$','Interpreter','latex');
ylabel('$$f(\psi)$$','Interpreter','latex');
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
grid on;
tuneFigure;
%saveFigure;