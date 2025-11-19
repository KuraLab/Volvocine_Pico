% phi の範囲を -pi から pi まで取得
phi = linspace(-pi, pi, 1000);

% 関数値
y = 0.5*sin(phi) - 1.0*sin(2*phi);

% プロット
figure;
plot(phi, y, 'LineWidth', 1.5);
xlabel('$$\phi$$','Interpreter','latex');
ylabel('$$f(\phi)$$','Interpreter','latex');
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
grid on;
tuneFigure;
saveFigure;