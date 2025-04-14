clear;
N = 5;
dt = 0.01;
time = 20;
steps = time/dt;
q = zeros(steps, N*2);
Order = zeros(steps, 1);
Theta = zeros(steps, 1);
Omega = 2*pi*ones(1, N);
kappa = 5*ones(1, N);
Sigmas = zeros(steps, N);
thetas = 0:2*pi/N:2*pi-2*pi/N;
shift = 0;
position = [100 100 500 400];

c1 = 1;      % 回転方向
c2 = 0.1;    % 平進方向

efficiency = ones(1, N); % 各モジュールの発揮率 (基本は1)
efficiency(3) = 0;       % 例: モジュール2が故障している場合

q(1, 1:N) = 2*pi*(rand(1, N)-0.5);

for i = 2:steps
    for j = 1:N
        Sigma = 0;
        for k = 1:N
            % calculateSigma関数を呼び出し
            Sigma = Sigma - (calculateSigma(k, j, N, c1, c2, thetas, efficiency) * ...
                sin(mod(q(i-1, j)-q(i-1, k)+pi, 2*pi)-pi));
        end
        Sigmas(i-1, j) = Sigma/N;
        q(i-1, j+N) = Omega(1, j) - kappa(1, j)*(Sigma/N);
        q(i, j) = q(i-1, j) + q(i-1, j+N)*dt;
        q(i, j) = mod(q(i, j)+pi, 2*pi)-pi;
    end
end

for i = 1:steps
    x = sum(sin(q(i, 1:N)))/N;
    y = sum(cos(q(i, 1:N)))/N;
    Order(i, 1) = sqrt(x^2 + y^2);
    Theta(i, 1) = mean(mod(q(i, 1)+pi, 2*pi)-pi);
end

% 発揮率が0でないモジュールのインデックスを取得
activeModules = find(efficiency > 0);

% 位相差の計算
phaseDiff = mod(q(:, activeModules) - Theta + pi + 0.05, 2*pi) - pi - shift;

% 値がpi以上ジャンプする箇所にNaNを挿入
for j = 1:length(activeModules)
    for i = 2:steps
        if abs(phaseDiff(i, j) - phaseDiff(i-1, j)) >= pi
            phaseDiff(i, j) = NaN;
        end
    end
end

% 位相差のプロット
figure('Position', position);
plot(dt:dt:steps*dt, phaseDiff)
grid on
set(gca, 'TickLabelInterpreter', 'latex')
xlim([0, steps*dt])
ylim([-pi+shift pi+shift])
ylabel('$$\phi_j - \phi_1$$');
yticks([-pi -pi/2 0 pi/2 pi 3*pi/2]);
yticklabels({'$$-\pi$$', '$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$', '$\pi$', '$3\pi/2$'})
xlabel('Time $$t$$ [s]');
tuneFigure;

% 秩序パラメータのプロット
figure('Position', position);
hold on
grid on
plot(dt:dt:steps*dt, Order)
ylim([0 1])
xlim([0, steps*dt])
ylabel('r');
xlabel('Time $$t$$ [s]');
tuneFigure;

% 秩序パラメータ計算 (OrderEXの処理を統合)
order = zeros(steps, 6);
P = zeros(2, N);
for i = 1:N
    P(1, i) = sin(thetas(i));
    P(2, i) = -cos(thetas(i));
end

for i = 1:steps
    for j = 1:N
        order(i, 1) = order(i, 1) + efficiency(j) * P(1, j)*sin(q(i, j))/N;
        order(i, 2) = order(i, 2) + efficiency(j) * P(1, j)*cos(q(i, j))/N;
        order(i, 4) = order(i, 4) + efficiency(j) * P(2, j)*sin(q(i, j))/N;
        order(i, 5) = order(i, 5) + efficiency(j) * P(2, j)*cos(q(i, j))/N;
    end
    order(i, 3) = sqrt(order(i, 1)^2 + order(i, 2)^2);
    order(i, 6) = sqrt(order(i, 4)^2 + order(i, 5)^2);
end

% 統合された秩序パラメータのプロット
figure;
plot(dt:dt:steps*dt, order(:, 3), dt:dt:steps*dt, order(:, 6))
ylim([0 1]);
xlim([0 steps*dt])
ylabel("R")
xlabel("Time t[s]")
legend("$$R_1$$", "$$R_2$$")
grid on
tuneFigure;
lineHandles = findobj(gca, 'Type', 'line');
set(lineHandles, 'LineWidth', 5);
set(gca, 'FontSize', 24);

% calculateSigma関数
function result = calculateSigma(k, j, N, c1, c2, thetas, efficiency)
    result = efficiency(k) * (c2 + c1 * cos(thetas(mod(k - j, N) + 1)));
end