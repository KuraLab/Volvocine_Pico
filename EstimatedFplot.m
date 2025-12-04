%% ============================
%  Spring 2 データ
%% ============================
freq2 = [230 235 240 245 250 255 260 265];
freq2_norm = (freq2 - 250) / 100;

phi2_a = [-0.790055, -0.666673, -0.453907, -0.283704, ...
           0.049815,  0.125068,  0.411129,  0.641414];

phi2_b = [-1.36336, -1.58014, -2.01999, -2.26026, ...
          -3.07159,  2.32331,  1.86263,  1.56455];

phi2_all  = [phi2_a, phi2_b];
freq2_all = [freq2,   freq2];
freq2_all_norm = (freq2_all - 250) / 100;


%% ============================
%  Spring 1 データ（手書きから読み取り）
%% ============================
% 周波数 240:   -0.51722,  -1.032
%        245:   -0.349878, 2.69335, 1.35462, -6.9891
%        250:   -0.195281, 3.0993,  1.18133, -0.902885
%        255:    0.015736,-2.91941, 0.976656,-0.793484
%        260:    0.266414, 0.463211

freq1_all = [...
    240 240, ...                         % 2 点
    245 245 245 245, ...                 % 4 点
    250 250 250 250, ...                 % 4 点
    255 255 255 255, ...                 % 4 点
    260 260];                            % 2 点

phi1_all  = [...
   -0.51722,  -1.032, ...
   -0.349878, 2.69335, 1.35462, -6.9891, ...
   -0.195281, 3.0993,  1.18133, -0.902885, ...
    0.015736,-2.91941, 0.976656,-0.793484, ...
    0.266414, 0.463211];
freq1_norm = (freq1_all - 250) / 100;

freq1_counts = [2, 4, 4, 4, 2];
freq1_stable = [];
freq1_unstable = [];
phi1_stable = [];
phi1_unstable = [];

cursor = 1;
for i = 1:numel(freq1_counts)
    n = freq1_counts(i);
    idx = cursor:(cursor + n - 1);
    split = n / 2;
    stable_idx = idx(1:split);
    unstable_idx = idx(split+1:end);

    freq1_stable = [freq1_stable, freq1_all(stable_idx)]; %#ok<AGROW>
    freq1_unstable = [freq1_unstable, freq1_all(unstable_idx)]; %#ok<AGROW>
    phi1_stable = [phi1_stable, phi1_all(stable_idx)]; %#ok<AGROW>
    phi1_unstable = [phi1_unstable, phi1_all(unstable_idx)]; %#ok<AGROW>

    cursor = cursor + n;
end

freq1_stable_norm = (freq1_stable - 250) / 100;
freq1_unstable_norm = (freq1_unstable - 250) / 100;


%% ============================
%  フィット関数定義
%  f(phi) = a*sin(phi+alpha) + b*sin(2phi+beta) + c
%% ============================

fitfun = @(params,phi) ...
    params(1).*sin(phi + params(2)) + ...
    params(3).*sin(2*phi + params(4)) + ...
    params(5);        % params = [a, alpha, b, beta, c]


%% ============================
%  Spring 2 のフィット
%% ============================

p0_2 = [10, 0, 10, 0, mean(freq2_all)];
lb   = [-Inf, -2*pi, -Inf, -2*pi, -Inf];
ub   = [ Inf,  2*pi,  Inf,  2*pi,  Inf];

opts = optimoptions('lsqcurvefit','Display','off');
params2 = lsqcurvefit(fitfun, p0_2, phi2_all, freq2_all, lb, ub, opts);


%% ============================
%  Spring 1 のフィット
%% ============================

p0_1 = [10, 0, 10, 0, mean(freq1_all)];
params1 = lsqcurvefit(fitfun, p0_1, phi1_all, freq1_all, lb, ub, opts);


%% ============================
%  フィット曲線生成
%% ============================

phi_fit  = linspace(-pi, pi, 500);
freq_fit1 = fitfun(params1, phi_fit);
freq_fit2 = fitfun(params2, phi_fit);
freq_fit1_norm = (freq_fit1 - 250) / 100;
freq_fit2_norm = (freq_fit2 - 250) / 100;


%% ============================
%  プロット（Spring1 & Spring2 を重ねる）
%% ============================


% ---- Spring 2 ----
figure;
scatter(phi2_a, freq2_norm, 70, 'o', 'filled', ...
    'MarkerFaceAlpha',0.8, 'DisplayName','Stable equilibrium'); hold on;
scatter(phi2_b, freq2_norm, 70, 's', 'filled', ...
    'MarkerFaceAlpha',0.8, 'DisplayName','Unstable equilibrium');
plot(phi_fit, freq_fit2_norm, 'LineWidth',2, 'DisplayName','Spring 2 fit');
apply_axis_format(gca, freq2_all_norm);
legend({'Stable equilibrium','Unstable equilibrium','Spring 2 fit'}, 'Location','best');
xlabel('$$\psi$$','Interpreter','latex');
ylabel('$$f(\psi)/\pi$$','Interpreter','latex');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','on');
tuneFigure;
saveFigure;

% ---- Spring 1 ----
figure;
scatter(phi1_stable, freq1_stable_norm, 70, 'o', 'filled', ...
    'MarkerFaceAlpha',0.8, 'DisplayName','Spring 1 stable'); hold on;
scatter(phi1_unstable, freq1_unstable_norm, 70, 's', 'filled', ...
    'MarkerFaceAlpha',0.8, 'DisplayName','Spring 1 unstable');
%plot(phi_fit, freq_fit1_norm, 'LineWidth',2,'DisplayName','Spring 1 fit');
apply_axis_format(gca, [freq1_stable_norm, freq1_unstable_norm]);
legend({'Spring 1 stable','Spring 1 unstable','Spring 1 fit'}, 'Location','best');
xlabel('$$\psi$$','Interpreter','latex');
ylabel('$$f(\psi)/\pi$$','Interpreter','latex');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','on');
tuneFigure;
saveFigure;

function apply_axis_format(ax, freq_vals)
    axes(ax); %#ok<LAXES>
    grid on;
    xlim([-pi, pi]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    fmin = min(freq_vals);
    fmax = max(freq_vals);
    pad  = max(1e-6, (fmax - fmin)*0.10);
    ylim([fmin - pad, fmax + pad]);
end
