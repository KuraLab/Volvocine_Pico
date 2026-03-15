freq_center = 250;
freq_scale = 100;
y_scale = pi;
clear_existing_figures = true;
save_figures = false;
show_fit_curves = false;
scale_freq = @(freq_vals) y_scale * (freq_vals - freq_center) / freq_scale;
scale_coeff = @(coeff_vals) y_scale * coeff_vals / freq_scale;
y_axis_label = '$$f(\psi)$$';

if clear_existing_figures
    close all hidden force;
    existing_figures = findall(groot, 'Type', 'figure');
    if ~isempty(existing_figures)
        delete(existing_figures);
    end
    drawnow;
end

%% ============================
%  Spring 2 データ
%% ============================
freq2 = [230 235 240 245 250 255 260 265];
freq2_norm = scale_freq(freq2);

phi2_a = [-0.790055, -0.666673, -0.453907, -0.283704, ...
           0.049815,  0.125068,  0.411129,  0.641414];

phi2_b = [-1.36336, -1.58014, -2.01999, -2.26026, ...
          -3.07159,  2.32331,  1.86263,  1.56455];

% 位相差の符号が反転していたため補正（全データの符号をひっくり返す）
phi2_a = -phi2_a;
phi2_b = -phi2_b;

phi2_all  = [phi2_a, phi2_b];
freq2_all = [freq2,   freq2];
freq2_all_norm = scale_freq(freq2_all);

%% ============================
%  Spring 3 データ
%% ============================
freq3 = [230 235 240 245 250 255 260 265];
freq3_norm = scale_freq(freq3);

phi3_a = [-0.861, -0.637, -0.475, -0.282, ...
           0.000,  0.246,  0.474,  0.669];

phi3_b = [-1.226, -2.133, -2.612, -2.787, ...
          -2.946,  2.254,  1.356,  1.114];

phi3_a = -phi3_a;
phi3_b = -phi3_b;

phi3_all  = [phi3_a, phi3_b];
freq3_all = [freq3,   freq3];
freq3_all_norm = scale_freq(freq3_all);

%% ============================
%  Spring 5 データ
%% ============================
freq5 = [235 240 245 250 255 260 265 270];
freq5_norm = scale_freq(freq5);

phi5_a = [-0.969, -0.687, -0.389, -0.090, ...
           0.439,  0.672,  0.871,  1.257];

phi5_b = [-1.934, -2.789, -3.038, 3.035, ...
          2.688,  2.226,  2.178,  1.957];

phi5_a = -phi5_a;
phi5_b = -phi5_b;

phi5_all  = [phi5_a, phi5_b];
freq5_all = [freq5,   freq5];
freq5_all_norm = scale_freq(freq5_all);

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
   -0.349878, 2.69335, 1.35462, -0.9891, ...
   -0.195281, 3.0993,  1.18133, -1.402885, ...
    0.015736,-2.91941, 0.976656,-1.793484, ...
    0.266414, 0.463211];
phi1_all = -phi1_all;
freq1_norm = scale_freq(freq1_all);

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

freq1_stable_norm = scale_freq(freq1_stable);
freq1_unstable_norm = scale_freq(freq1_unstable);


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
%  Spring 3 のフィット
%% ============================

p0_3 = [10, 0, 10, 0, mean(freq3_all)];
params3 = lsqcurvefit(fitfun, p0_3, phi3_all, freq3_all, lb, ub, opts);


%% ============================
%  Spring 5 のフィット
%% ============================

p0_5 = [10, 0, 10, 0, mean(freq5_all)];
params5 = lsqcurvefit(fitfun, p0_5, phi5_all, freq5_all, lb, ub, opts);


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
freq_fit3 = fitfun(params3, phi_fit);
freq_fit5 = fitfun(params5, phi_fit);
freq_fit1_norm = scale_freq(freq_fit1);
freq_fit2_norm = scale_freq(freq_fit2);
freq_fit3_norm = scale_freq(freq_fit3);
freq_fit5_norm = scale_freq(freq_fit5);


%% ============================
%  プロット（Spring1 & Spring2 を重ねる）
%% ============================

% Consistent colors per spring (one color per spring)
c1 = [0.1216, 0.4667, 0.7059]; % Spring 1
c2 = [1.0000, 0.4980, 0.0549]; % Spring 2
c3 = [0.1725, 0.6275, 0.1725]; % Spring 3
c5 = [0.8392, 0.1529, 0.1569]; % Spring 5


% ---- Spring 2 ----
figure;
scatter(phi2_a, freq2_norm, 70, 'o', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c2, 'MarkerEdgeColor',c2, 'DisplayName','Stable equilibrium'); hold on;
scatter(phi2_b, freq2_norm, 70, 's', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c2, 'MarkerEdgeColor',c2, 'DisplayName','Unstable equilibrium');
maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit2_norm, 'LineWidth',2, 'Color',c2, 'DisplayName','Spring 2 fit');
apply_axis_format(gca, freq2_all_norm);
legend_labels = {'Stable equilibrium','Unstable equilibrium'};
if show_fit_curves
    legend_labels{end+1} = 'Spring 2 fit'; %#ok<AGROW>
end
legend(legend_labels, 'Location','best');
xlabel('$$\psi$$','Interpreter','latex');
ylabel(y_axis_label,'Interpreter','latex');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','on');
tuneFigure;
maybe_save_figure(save_figures);

%% ============================
%  Compare fit coefficients (1st vs 2nd harmonic)
%% ============================

springs = {'Spring 1','Spring 2','Spring 3','Spring 5'};
a1 = [params1(1), params2(1), params3(1), params5(1)]; % coefficient of sin(psi + alpha)
b2 = [params1(3), params2(3), params3(3), params5(3)]; % coefficient of sin(2psi + beta)
a1_plot = scale_coeff(a1);
b2_plot = scale_coeff(b2);

T = table(springs(:), a1(:), b2(:), abs(b2(:))./max(1e-12,abs(a1(:))), ...
    'VariableNames', {'Spring','a1_sin1','b2_sin2','abs_b2_over_abs_a1'});
disp('=== Fit coefficient comparison ===');
disp(T);

figure;
bar([a1_plot(:), b2_plot(:)]);
grid on;
set(gca, 'XTick', 1:numel(springs), 'XTickLabel', springs, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
xlabel('Spring');
ylabel('$$f(\psi)$$','Interpreter','latex');
legend({'a (sin(\psi+\alpha))','b (sin(2\psi+\beta))'}, 'Location', 'best');
title('Fit coefficients for 1st and 2nd harmonic terms');
tuneFigure;
maybe_save_figure(save_figures);

% ---- Spring 3 ----
figure;
scatter(phi3_a, freq3_norm, 70, 'o', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c3, 'MarkerEdgeColor',c3, 'DisplayName','Stable equilibrium'); hold on;
scatter(phi3_b, freq3_norm, 70, 's', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c3, 'MarkerEdgeColor',c3, 'DisplayName','Unstable equilibrium');
maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit3_norm, 'LineWidth',2, 'Color',c3, 'DisplayName','Spring 3 fit');
apply_axis_format(gca, freq3_all_norm);
legend_labels = {'Stable equilibrium','Unstable equilibrium'};
if show_fit_curves
    legend_labels{end+1} = 'Spring 3 fit'; %#ok<AGROW>
end
legend(legend_labels, 'Location','best');
xlabel('$$\psi$$','Interpreter','latex');
ylabel(y_axis_label,'Interpreter','latex');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','on');
tuneFigure;
maybe_save_figure(save_figures);

% ---- Spring 5 ----
figure;
scatter(phi5_a, freq5_norm, 70, 'o', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c5, 'MarkerEdgeColor',c5, 'DisplayName','Stable equilibrium'); hold on;
scatter(phi5_b, freq5_norm, 70, 's', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c5, 'MarkerEdgeColor',c5, 'DisplayName','Unstable equilibrium');
maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit5_norm, 'LineWidth',2, 'Color',c5, 'DisplayName','Rigid fit');
apply_axis_format(gca, freq5_all_norm);
legend_labels = {'Stable equilibrium','Unstable equilibrium'};
if show_fit_curves
    legend_labels{end+1} = 'Rigid fit'; %#ok<AGROW>
end
legend(legend_labels, 'Location','best');
xlabel('$$\psi$$','Interpreter','latex');
ylabel(y_axis_label,'Interpreter','latex');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','on');
tuneFigure;
maybe_save_figure(save_figures);

% ---- Spring 1 ----
figure;
scatter(phi1_stable, freq1_stable_norm, 70, 'o', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c1, 'MarkerEdgeColor',c1, 'DisplayName','Stable equilibrium'); hold on;
scatter(phi1_unstable, freq1_unstable_norm, 70, 's', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerFaceColor',c1, 'MarkerEdgeColor',c1, 'DisplayName','Unstable equilibrium');
maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit1_norm, 'LineWidth',2, 'Color',c1, 'DisplayName','Spring 1 fit');
apply_axis_format(gca, [freq1_stable_norm, freq1_unstable_norm]);
legend_labels = {'Spring 1 stable','Spring 1 unstable'};
if show_fit_curves
    legend_labels{end+1} = 'Spring 1 fit'; %#ok<AGROW>
end
legend(legend_labels, 'Location','best');
xlabel('$$\psi$$','Interpreter','latex');
ylabel(y_axis_label,'Interpreter','latex');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','on');
tuneFigure;
maybe_save_figure(save_figures);

%% ============================
%  Overlay plot (Spring 1/2/3/5)
%% ============================

figure;

% Spring 1
scatter(phi1_stable, freq1_stable_norm, 60, 'o', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c1, 'MarkerEdgeColor',c1, 'DisplayName','Stable fixed point'); hold on;
scatter(phi1_unstable, freq1_unstable_norm, 60, 's', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c1, 'MarkerEdgeColor',c1, 'DisplayName','Unstable fixed point');
h1 = maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit1_norm, 'LineWidth',2, 'Color',c1, 'DisplayName','Spring 1 fit');

% Spring 2
scatter(phi2_a, freq2_norm, 60, 'o', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c2, 'MarkerEdgeColor',c2, 'DisplayName','Stable fixed point');
scatter(phi2_b, freq2_norm, 60, 's', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c2, 'MarkerEdgeColor',c2, 'DisplayName','Unstable fixed point');
h2 = maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit2_norm, 'LineWidth',2, 'Color',c2, 'DisplayName','Spring 2 fit');

% Spring 3
scatter(phi3_a, freq3_norm, 60, 'o', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c3, 'MarkerEdgeColor',c3, 'DisplayName','Stable fixed point');
scatter(phi3_b, freq3_norm, 60, 's', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c3, 'MarkerEdgeColor',c3, 'DisplayName','Unstable fixed point');
h3 = maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit3_norm, 'LineWidth',2, 'Color',c3, 'DisplayName','Spring 3 fit');

% Spring 5
scatter(phi5_a, freq5_norm, 60, 'o', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c5, 'MarkerEdgeColor',c5, 'DisplayName','Stable fixed point');
scatter(phi5_b, freq5_norm, 60, 's', 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerFaceColor',c5, 'MarkerEdgeColor',c5, 'DisplayName','Unstable fixed point');
h5 = maybe_plot_fit_curve(show_fit_curves, phi_fit, freq_fit5_norm, 'LineWidth',2, 'Color',c5, 'DisplayName','Spring 5 fit');

apply_axis_format(gca, [freq1_stable_norm, freq1_unstable_norm, freq2_all_norm, freq3_all_norm, freq5_all_norm]);
if show_fit_curves
    legend([h1, h2, h3, h5], {'Spring 1 fit','Spring 2 fit','Spring 3 fit','Spring 5 fit'}, 'Location','best');
end
xlabel('$$\psi$$','Interpreter','latex');
ylabel(y_axis_label,'Interpreter','latex');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','on');
tuneFigure;
maybe_save_figure(save_figures);

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

function maybe_save_figure(save_figures)
    if save_figures
        saveFigure;
    end
end

function plot_handle = maybe_plot_fit_curve(show_fit_curves, x_vals, y_vals, varargin)
    if show_fit_curves
        plot_handle = plot(x_vals, y_vals, varargin{:});
    else
        plot_handle = gobjects(0);
    end
end
