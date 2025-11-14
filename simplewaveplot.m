% phi の範囲
phi = linspace(0, 2*pi, 1000);

% 関数値
y = -sin(phi) - 1.5*sin(2*phi);

% プロット
figure;
plot(phi, y, 'LineWidth', 1.5);
xlabel('\phi');
ylabel('sin(\phi) + sin(2\phi)');
title('Plot of sin(\phi) + sin(2\phi)');
grid on;
